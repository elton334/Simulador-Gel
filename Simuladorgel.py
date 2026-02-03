import streamlit as st
import plotly.graph_objects as go
import math
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO, BytesIO
from Bio import SeqIO

# --- CONFIGURAÇÃO INICIAL ---
st.set_page_config(page_title="Simulador de Gel Interativo", layout="wide", page_icon="🧬")

# Converte enzimas para string
TODAS_ENZIMAS = sorted([str(e) for e in CommOnly])

# Dados de Ladders
LADDERS = {
    "1kb Plus DNA Ladder": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "1kb DNA Ladder (Genérico)": [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
    "100bp DNA Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}

def processar_upload(input_data):
    """Lê arquivos FASTA (.fasta, .txt) ou SnapGene (.dna)."""
    try:
        nome_arquivo = input_data.name
        # Remove extensão para usar como nome sugerido
        nome_sugerido = nome_arquivo.rsplit('.', 1)[0]
        
        # --- CASO 1: Arquivo SnapGene (.dna) ---
        if nome_arquivo.lower().endswith('.dna'):
            try:
                bytes_io = BytesIO(input_data.getvalue())
                record = SeqIO.read(bytes_io, "snapgene")
                return nome_sugerido, str(record.seq).upper()
            except Exception as e:
                return "Erro", f"Erro .dna: {str(e)}"

        # --- CASO 2: Arquivo Texto (FASTA/TXT) ---
        bytes_data = input_data.getvalue()
        try:
            conteudo = bytes_data.decode("utf-8")
        except UnicodeDecodeError:
            conteudo = bytes_data.decode("latin-1")

        if ">" in conteudo:
            try:
                iterator = SeqIO.parse(StringIO(conteudo), "fasta")
                record = next(iterator)
                # Usa ID do fasta ou nome do arquivo
                return record.id if record.id else nome_sugerido, str(record.seq).upper()
            except:
                pass 

        # Fallback: Texto Cru (Raw)
        linhas = conteudo.splitlines()
        seq_limpa = ""
        for linha in linhas:
            linha = linha.strip()
            if not linha or linha.startswith(">") or linha.startswith(";"): continue
            seq_limpa += linha
        
        seq_final = "".join(seq_limpa.split()).upper()
        if any(c not in "ATGCNRYKMSWBDHV" for c in seq_final[:100]): 
             return "Erro", "Arquivo inválido."

        return nome_sugerido, seq_final

    except Exception as e:
        return "Erro", str(e)

def processar_texto_manual(texto):
    try:
        if ">" in texto:
            iterator = SeqIO.parse(StringIO(texto), "fasta")
            record = next(iterator)
            return record.id, str(record.seq).upper()
        else:
            return "Seq Manual", "".join(texto.split()).upper()
    except:
        return "Erro", ""

def calcular_digestao(sequencia, enzimas, eh_circular):
    if not sequencia or sequencia.startswith("Erro"): return []
    sequencia = "".join([c for c in sequencia if c in "ATGCMRWSYKVHDBN"])
    if not sequencia: return []

    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)
    
    if eh_circular and not enzimas:
        return [
            (tamanho_total * 1.4, "Nicked", tamanho_total),
            (tamanho_total * 0.7, "Supercoiled", tamanho_total)
        ]
    
    if not enzimas: 
        return [(tamanho_total, "Linear", tamanho_total)]
    
    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))
    
    if not locais: 
        tipo = "Circ. (Sítio Ausente)" if eh_circular else "Lin. (Não Cortado)"
        return [(tamanho_total, tipo, tamanho_total)]
        
    fragmentos = []
    if not eh_circular:
        prev = 0
        for cut in locais:
            fragmentos.append(cut - prev)
            prev = cut
        fragmentos.append(tamanho_total - prev)
    else:
        if len(locais) == 1:
            fragmentos.append(tamanho_total)
        else:
            for i in range(len(locais)-1):
                fragmentos.append(locais[i+1] - locais[i])
            fragmentos.append((tamanho_total - locais[-1]) + locais[0])
            
    return [(frag, "Fragmento", frag) for frag in sorted(fragmentos, reverse=True)]

# --- INTERFACE ---
st.title("🧪 Simulador de Eletroforese Interativo")

with st.sidebar:
    st.header("Configurações")
    num_pocos = st.slider("Número de Poços", 1, 15, 3) 
    st.divider()
    inverter_cores = st.toggle("Inverter Cores", value=False)
    
    # SUGESTÃO DE OURO: Ajuste da % de Agarose (Muda o Zoom Vertical)
    agarose = st.slider("Concentração de Agarose (%)", 0.5, 2.0, 1.0, 0.1)
    st.caption("Ajustar a agarose altera a faixa de visualização (Zoom vertical).")

dados_para_plotar = []
labels_eixo_x = []
nomes_ladders = [] 
detalhes_hover = [] 

cols = st.columns(2)

for i in range(num_pocos):
    col_atual = cols[i % 2]
    with col_atual:
        with st.expander(f"Poço {i+1}", expanded=(i==0)):
            tipo = st.radio(f"Conteúdo {i+1}:", ["Amostra", "Ladder"], key=f"t_{i}", horizontal=True)
            
            # --- SUGESTÃO IMPLEMENTADA: NOME DO POÇO ---
            # O padrão é o número, mas você pode mudar para "pUC19", "Clone 1", etc.
            rotulo_padrao = str(i+1)
            
            if tipo == "Ladder":
                lad = st.selectbox("Ladder:", list(LADDERS.keys()), key=f"l_{i}")
                ladder_data = [(tam, "Ladder", tam) for tam in LADDERS[lad]]
                dados_para_plotar.append(ladder_data)
                
                # Se for ladder, sugerimos "M" ou o nome do ladder
                rotulo_custom = st.text_input("Nome no Gel:", value="M", key=f"lbl_{i}")
                labels_eixo_x.append(rotulo_custom)
                
                nomes_ladders.append(lad)
                detalhes_hover.append(lad)
            else:
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs(["Arquivo", "Texto"])
                seq, nome_arquivo = "", ""
                
                with tab_f:
                    up = st.file_uploader("Arquivo", type=['dna', 'fasta', 'txt', 'fa'], key=f"u_{i}")
                    if up: 
                        nome_arquivo, seq = processar_upload(up)
                        if nome_arquivo == "Erro": 
                            st.error(seq); seq = ""
                with tab_t:
                    txt = st.text_area("Colar Sequência", height=70, key=f"tx_{i}")
                    if txt and not seq: 
                        nome_t, seq_t = processar_texto_manual(txt)
                        if nome_t != "Seq Manual": nome_arquivo = nome_t
                        seq = seq_t
                
                c1, c2 = st.columns(2)
                circ = c1.checkbox("Circular?", True, key=f"c_{i}")
                enz = c2.multiselect("Enzimas", TODAS_ENZIMAS, key=f"e_{i}")
                
                # Usa o nome do arquivo como sugestão de rótulo se disponível
                val_rotulo = nome_arquivo if nome_arquivo else str(i+1)
                rotulo_custom = st.text_input("Nome no Gel:", value=val_rotulo[:10], key=f"lbl_{i}")
                labels_eixo_x.append(rotulo_custom)

                info_texto = f"Circular: {circ}<br>Enzimas: {', '.join(enz) if enz else 'Uncut'}"
                detalhes_hover.append(info_texto)

                if seq:
                    try:
                        res = calcular_digestao(seq, enz, circ)
                        dados_para_plotar.append(res)
                    except Exception as e:
                        dados_para_plotar.append([])
                else:
                    dados_para_plotar.append([])

st.divider()

if any(dados_para_plotar):
    bg_color = 'white' if inverter_cores else '#1e1e1e'
    line_color = 'black' if inverter_cores else 'white'
    text_color = 'black' if inverter_cores else 'white'
    
    # Define limites verticais baseado na agarose (Física aproximada)
    # 0.7% vê até 20kb. 2.0% vê até 2kb.
    min_view = 50 + (100 * (agarose - 0.5)) 
    max_view = 25000 / (agarose * 0.8)

    fig = go.Figure()

    for i, lista_bandas in enumerate(dados_para_plotar):
        x_center = i + 1
        eh_ladder = (nomes_ladders[i] is not None)
        
        if lista_bandas:
             massa_total = sum([b[2] for b in lista_bandas]) if not eh_ladder else 1
        
        for (tam_aparente, tipo_banda, tam_real) in lista_bandas:
            # Filtra bandas fora do alcance visual da agarose
            if tam_aparente < min_view or tam_aparente > max_view: continue

            width = 2; opacity = 0.8
            if eh_ladder:
                if tam_aparente in [3000, 1000, 500]: width = 7; opacity = 1.0
                elif tam_aparente >= 5000: width = 5; opacity = 0.9
                else: width = 3; opacity = 0.7
            else:
                if tipo_banda == "Supercoiled": fracao = 0.7
                elif tipo_banda == "Nicked": fracao = 0.3
                else: fracao = tam_real / massa_total if massa_total > 0 else 0.5
                width = 3 + (8 * fracao)
                opacity = 0.6 + (0.4 * fracao)

            largura_banda = 0.28 
            
            fig.add_trace(go.Scatter(
                x=[x_center - largura_banda, x_center + largura_banda],
                y=[tam_aparente, tam_aparente],
                mode='lines+markers',
                line=dict(color=line_color, width=width),
                marker=dict(color=line_color, size=width, symbol='circle'),
                opacity=opacity,
                showlegend=False,
                hoverinfo='text',
                hovertext=f"<b>~{int(tam_aparente)} pb</b><br>Real: {tam_real}<br>{labels_eixo_x[i]}"
            ))

            if eh_ladder:
                fig.add_trace(go.Scatter(
                    x=[x_center - 0.45], y=[tam_aparente], mode="text",
                    text=[str(tam_aparente)], textposition="middle left",
                    textfont=dict(color=text_color, size=10),
                    showlegend=False, hoverinfo='skip'
                ))

    # --- LAYOUT FINAL ---
    LARGURA_MINIMA = 12
    max_range = max(num_pocos, LARGURA_MINIMA) + 0.5

    fig.update_layout(
        plot_bgcolor=bg_color, paper_bgcolor=bg_color,
        height=700, margin=dict(t=40, b=40, l=40, r=40),
        xaxis=dict(
            tickmode='array', tickvals=list(range(1, num_pocos + 1)),
            ticktext=labels_eixo_x, # AGORA USA SEUS NOMES PERSONALIZADOS
            tickfont=dict(color=text_color, size=14, family='Arial Black'),
            showgrid=False, zeroline=False, range=[0.2, max_range] 
        ),
        yaxis=dict(
            type='log',
            range=[math.log10(min_view), math.log10(max_view)], # Zoom dinâmico pela agarose
            showgrid=False, zeroline=False, showticklabels=False
        )
    )
    
    fig.add_annotation(x=-0.05, y=1, xref="paper", yref="paper", text="pb", showarrow=False, font=dict(color=text_color, size=14, family="Arial Black"))
    st.plotly_chart(fig, use_container_width=True)

else:
    st.info("Adicione amostras para gerar o gel.")
