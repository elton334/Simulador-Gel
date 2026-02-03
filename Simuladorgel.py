import streamlit as st
import plotly.graph_objects as go
import math
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO
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

def processar_fasta(input_data, is_file=False):
    try:
        if is_file:
            stringio = StringIO(input_data.getvalue().decode("utf-8"))
            iterator = SeqIO.parse(stringio, "fasta")
        else:
            if ">" in input_data:
                from io import StringIO
                iterator = SeqIO.parse(StringIO(input_data), "fasta")
            else:
                return "Seq Manual", "".join(input_data.split()).upper()
        
        record = next(iterator)
        return record.id, str(record.seq).upper()
    except Exception:
        return "Erro", ""

def calcular_digestao(sequencia, enzimas, eh_circular):
    if not sequencia: return []
    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)
    
    # CASO ESPECIAL: Plasmídeo Circular Não Cortado (Uncut)
    if eh_circular and not enzimas:
        # Simulação de migração anômala
        # Supercoiled: Corre mais rápido (aprox 0.7x do tamanho linear)
        # Nicked/Relaxed: Corre mais devagar (aprox 1.4x do tamanho linear)
        
        # Retorna lista de tuplas: (Tamanho Aparente, Tipo, Tamanho Real)
        return [
            (tamanho_total * 1.4, "Nicked (Relaxed)", tamanho_total),
            (tamanho_total * 0.7, "Supercoiled", tamanho_total)
        ]
    
    if not enzimas: 
        # Linear não cortado
        return [(tamanho_total, "Linear", tamanho_total)]
    
    # Digestão Normal
    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))
    
    if not locais: 
        # Enzima não cortou (site ausente)
        return [(tamanho_total, "Circular (Não Cortado - Sítio Ausente)", tamanho_total)] if eh_circular else [(tamanho_total, "Linear (Não Cortado)", tamanho_total)]
        
    fragmentos = []
    if not eh_circular:
        prev = 0
        for cut in locais:
            fragmentos.append(cut - prev)
            prev = cut
        fragmentos.append(tamanho_total - prev)
    else:
        if len(locais) == 1:
            fragmentos.append(tamanho_total) # Linearizou
        else:
            for i in range(len(locais)-1):
                fragmentos.append(locais[i+1] - locais[i])
            fragmentos.append((tamanho_total - locais[-1]) + locais[0])
            
    # Formata saida padrão para digestão: (Tamanho, Tipo, Real)
    return [(frag, "Fragmento", frag) for frag in sorted(fragmentos, reverse=True)]

# --- INTERFACE ---
st.title("🧪 Simulador de Eletroforese Interativo")
st.markdown("Passe o mouse sobre as bandas para ver detalhes.")

with st.sidebar:
    st.header("Configurações")
    num_pocos = st.slider("Número de Poços", 1, 15, 10)
    st.divider()
    inverter_cores = st.toggle("Inverter Cores (Modo Impressão)", value=False)
    st.caption("Nota: Plasmídeos não cortados mostram bandas Supercoiled e Nicked.")

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
            
            if tipo == "Ladder":
                lad = st.selectbox("Ladder:", list(LADDERS.keys()), key=f"l_{i}")
                # Padroniza ladder para o formato (Tamanho, Tipo, Real)
                ladder_data = [(tam, "Ladder", tam) for tam in LADDERS[lad]]
                dados_para_plotar.append(ladder_data)
                labels_eixo_x.append("M")
                nomes_ladders.append(lad)
                detalhes_hover.append(lad)
            else:
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs(["Arquivo", "Texto"])
                seq, nome = "", f"{i+1}"
                
                with tab_f:
                    up = st.file_uploader("FASTA", key=f"u_{i}")
                    if up: nome, seq = processar_fasta(up, True)
                with tab_t:
                    txt = st.text_area("Seq", height=70, key=f"tx_{i}")
                    if txt and not seq: 
                        nome_t, seq_t = processar_fasta(txt, False)
                        if nome_t != "Seq Manual": nome = nome_t
                        seq = seq_t
                
                c1, c2 = st.columns(2)
                circ = c1.checkbox("Circular?", True, key=f"c_{i}")
                enz = c2.multiselect("Enzimas", TODAS_ENZIMAS, key=f"e_{i}")
                
                info_texto = f"Circular: {circ}<br>Enzimas: {', '.join(enz) if enz else 'Nenhuma (Uncut)'}"
                detalhes_hover.append(info_texto)

                if seq:
                    try:
                        # Agora retorna tuplas (Tamanho Aparente, Tipo, Tamanho Real)
                        res = calcular_digestao(seq, enz, circ)
                        dados_para_plotar.append(res)
                        labels_eixo_x.append(str(i+1))
                    except:
                        dados_para_plotar.append([])
                        labels_eixo_x.append("Erro")
                else:
                    dados_para_plotar.append([])
                    labels_eixo_x.append(str(i+1))

st.divider()

if any(dados_para_plotar):
    bg_color = 'white' if inverter_cores else '#1e1e1e'
    line_color = 'black' if inverter_cores else 'white'
    text_color = 'black' if inverter_cores else 'white'
    
    fig = go.Figure()

    for i, lista_bandas in enumerate(dados_para_plotar):
        x_center = i + 1
        eh_ladder = (labels_eixo_x[i] == "M")
        
        # Calcula massa total SOMENTE dos fragmentos reais para proporção
        # Se for uncut (plasmídeo inteiro), a massa total é o tamanho do plasmídeo (aparece 2x mas a massa é dividida)
        # Para simplificar visualização: usamos o tamanho real da maior banda como referência de massa
        if lista_bandas:
             massa_total = sum([b[2] for b in lista_bandas]) if not eh_ladder else 1
        
        for (tam_aparente, tipo_banda, tam_real) in lista_bandas:
            # --- ESTÉTICA ---
            width = 2
            opacity = 0.8
            
            if eh_ladder:
                if tam_aparente in [3000, 1000, 500]: 
                    width = 7; opacity = 1.0
                elif tam_aparente >= 5000:
                    width = 5; opacity = 0.9
                else:
                    width = 3; opacity = 0.7
            else:
                # Se for plasmídeo uncut, as bandas Supercoiled e Nicked dividem a massa
                # Mas visualmente, Supercoiled costuma ser mais forte e nítida.
                if tipo_banda == "Supercoiled":
                    fracao = 0.7 # Supercoiled brilha mais pois é compacta
                elif tipo_banda == "Nicked (Relaxed)":
                    fracao = 0.3 # Nicked é mais difusa
                else:
                    fracao = tam_real / massa_total if massa_total > 0 else 0.5
                
                width = 3 + (8 * fracao)
                opacity = 0.6 + (0.4 * fracao)

            # --- DESENHO DA BANDA ---
            fig.add_trace(go.Scatter(
                x=[x_center - 0.35, x_center + 0.35],
                y=[tam_aparente, tam_aparente],
                mode='lines+markers',
                line=dict(color=line_color, width=width),
                marker=dict(color=line_color, size=width, symbol='circle'),
                opacity=opacity,
                showlegend=False,
                hoverinfo='text',
                # Tooltip Educativo: Mostra o tamanho aparente E o real
                hovertext=f"<b>Aparente: ~{int(tam_aparente)} pb</b><br>Real: {tam_real} pb<br>Tipo: {tipo_banda}<br>Poço: {labels_eixo_x[i]}"
            ))

            # --- RÓTULOS DO LADDER ---
            if eh_ladder:
                fig.add_trace(go.Scatter(
                    x=[x_center - 0.5], 
                    y=[tam_aparente],
                    mode="text",
                    text=[str(tam_aparente)],
                    textposition="middle left",
                    textfont=dict(color=text_color, size=10),
                    showlegend=False,
                    hoverinfo='skip'
                ))

    # --- LAYOUT ---
    fig.update_layout(
        plot_bgcolor=bg_color,
        paper_bgcolor=bg_color,
        height=700,
        margin=dict(t=40, b=40, l=40, r=40),
        xaxis=dict(
            tickmode='array',
            tickvals=list(range(1, num_pocos + 1)),
            ticktext=labels_eixo_x,
            tickfont=dict(color=text_color, size=14, family='Arial Black'),
            showgrid=False, zeroline=False, range=[0.2, num_pocos + 0.8] 
        ),
        yaxis=dict(
            type='log',
            range=[math.log10(80), math.log10(20000)], # Invertido
            showgrid=False, zeroline=False, showticklabels=False
        )
    )
    
    fig.add_annotation(x=-0.05, y=1, xref="paper", yref="paper", text="pb", showarrow=False, font=dict(color=text_color, size=14, family="Arial Black"))
    st.plotly_chart(fig, use_container_width=True)

else:
    st.info("Adicione amostras para gerar o gel.")
