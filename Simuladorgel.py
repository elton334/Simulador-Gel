import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math
import re
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO, BytesIO
from Bio import SeqIO

# --- 1. CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(
    page_title="BioSpark Studio Pro",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

# --- 2. GERENCIAMENTO DE ESTADO (O CÉREBRO DA v2.0) ---
# Inicializa variáveis para persistência se não existirem
if 'lab_notebook' not in st.session_state:
    st.session_state.lab_notebook = ""

# Função para injetar dados do Upload em Lote nos poços individuais
def handle_batch_upload():
    uploaded_files = st.session_state.batch_uploader
    if uploaded_files:
        for idx, file in enumerate(uploaded_files):
            # Limite de 15 poços
            if idx >= 15: break
            
            # Processa o arquivo usando a função existente
            name, seq = processar_upload(file)
            
            if name != "Erro":
                # Injeta diretamente no estado dos widgets dos poços
                # As chaves dos widgets são: tx_{i} para texto e lbl_{i} para rótulo
                st.session_state[f"tx_{idx}"] = seq
                st.session_state[f"lbl_{idx}"] = name[:15]
                # Define o tipo como 'Amostra' (Digestão) por padrão para facilitar
                # Nota: Mudar o radio button via estado é complexo, mantemos o padrão do usuário

# --- 3. SISTEMA DE TRADUÇÃO ---
TEXTS = {
    "header_title": { "PT": "BioSpark Studio", "EN": "BioSpark Studio" },
    "header_sub": { "PT": "Simulação Avançada: Digestão, PCR e Batch Processing.", "EN": "Advanced Simulation: Digestion, PCR and Batch Processing." },
    "sidebar_config": { "PT": "CONFIGURAÇÕES GERAIS", "EN": "GENERAL SETTINGS" },
    "sidebar_wells": { "PT": "Nº de Poços", "EN": "Well Count" },
    "sidebar_agarose": { "PT": "Agarose (%)", "EN": "Agarose (%)" },
    "batch_label": { "PT": "⚡ Upload em Lote (Arrastar múltiplos arquivos)", "EN": "⚡ Batch Upload (Drag multiple files)" },
    "notebook_title": { "PT": "📔 Caderno de Lab", "EN": "📔 Lab Notebook" },
    
    # Mantendo os textos clássicos
    "well_title": { "PT": "Poço", "EN": "Well" },
    "opt_sample": { "PT": "Digestão", "EN": "Digestion" },
    "opt_ladder": { "PT": "Ladder", "EN": "Ladder" },
    "opt_pcr": { "PT": "PCR", "EN": "PCR" },
    "sel_ladder": { "PT": "Selecione o Ladder:", "EN": "Select Ladder:" },
    "label_gel": { "PT": "Rótulo:", "EN": "Label:" },
    "tab_file": { "PT": "📂 Arquivo", "EN": "📂 File" },
    "tab_text": { "PT": "📝 Editor", "EN": "📝 Editor" },
    "upload_label": { "PT": "Upload único", "EN": "Single Upload" },
    "paste_label": { "PT": "Sequência", "EN": "Sequence" },
    "check_circular": { "PT": "Circular?", "EN": "Circular?" },
    "sel_enzymes": { "PT": "Enzimas", "EN": "Enzymes" },
    "pcr_fwd": { "PT": "Fwd (5'->3')", "EN": "Fwd (5'->3')" },
    "pcr_rev": { "PT": "Rev (3'->5' ou 5'->3')", "EN": "Rev (3'->5' or 5'->3')" },
    "result_title": { "PT": "Resultado da Eletroforese", "EN": "Electrophoresis Result" },
    "export_expander": { "PT": "Exportar Dados", "EN": "Export Data" },
    "btn_download": { "PT": "Baixar .csv", "EN": "Download .csv" },
    "empty_msg": { "PT": "Adicione amostras ou use o Upload em Lote.", "EN": "Add samples or use Batch Upload." },
    "created_by": { "PT": "Dev.", "EN": "Dev." },
    "report_bug": { "PT": "✉️ Reportar Erro", "EN": "✉️ Report Bug" },
    "warn_multiple": { "PT": "⚠️ Inespecífico!", "EN": "⚠️ Non-specific!" },
    "warn_no_product": { "PT": "Sem produto.", "EN": "No product." },
    "diag_fwd_fail": { "PT": "❌ Fwd fail", "EN": "❌ Fwd fail" },
    "diag_rev_fail": { "PT": "❌ Rev fail", "EN": "❌ Rev fail" },
    "ack_title": { "PT": "Apoio e Afiliação", "EN": "Support & Affiliation" }
}

# --- 4. ESTILO CSS (TURQUESA PRO) ---
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    .stApp { background: linear-gradient(180deg, #F0F9FF 0%, #FFFFFF 100%); font-family: 'Inter', sans-serif; }
    
    /* Sidebar */
    section[data-testid="stSidebar"] { background-color: #E0F7FA; border-right: 1px solid #B2EBF2; }
    section[data-testid="stSidebar"] .block-container { padding-top: 1rem; gap: 0.5rem; }

    /* Componentes */
    div[data-baseweb="slider"] div[class*="StyledThumb"] { background-color: #0F766E !important; border-color: #0F766E !important; }
    div[data-baseweb="slider"] div[class*="StyledTrack"] > div { background-color: #0F766E !important; }
    div[data-baseweb="slider"] div[class*="StyledTrack"] { background-color: #B2EBF2 !important; }
    
    h1, h2, h3 { color: #0F172A !important; font-weight: 700 !important; letter-spacing: -0.02em; }
    
    .stExpander { background-color: #FFFFFF; border-radius: 6px !important; border: 1px solid #E2E8F0 !important; box-shadow: 0 1px 2px 0 rgba(0,0,0,0.05) !important; margin-bottom: 0.5rem; }
    .stExpander:hover { border-color: #0F766E !important; }
    
    .stButton > button { background-color: #0F766E; color: white; border-radius: 6px; font-weight: 500; border: none; }
    .stButton > button:hover { background-color: #0d6e66; color: white; }
    
    span[data-baseweb="checkbox"] div { background-color: #0F766E !important; }
    button[data-baseweb="tab"] { font-size: 13px !important; padding: 5px 10px !important; }

    /* Área de Batch Upload */
    div[data-testid="stFileUploader"] { margin-bottom: 10px; }

    /* Rodapés */
    .sidebar-footer { margin-top: 15px; padding-top: 10px; border-top: 1px solid #CBD5E1; font-size: 11px; color: #333333; line-height: 1.4; }
    .sidebar-footer strong { color: #111827; font-weight: 600; }
    .bug-report { font-size: 10px; color: #4B5563; text-decoration: none; margin-left: 5px; }
    .bug-report:hover { color: #0F766E; text-decoration: underline; }
    .footer { width: 100%; text-align: center; padding: 20px 0; font-size: 11px; color: #64748B; border-top: 1px solid #CBD5E1; margin-top: 40px; opacity: 0.8; }
    
    .warning-text { color: #DC2626; font-weight: bold; font-size: 11px; margin: 0; }
    .error-text { color: #B91C1C; font-size: 11px; margin: 0; }
</style>
""", unsafe_allow_html=True)

# --- 5. LÓGICA BIOLÓGICA (MANTIDA V3.0) ---

TODAS_ENZIMAS = sorted([str(e) for e in CommOnly])
LADDERS = {
    "1kb Plus DNA Ladder": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "1kb DNA Ladder (Genérico)": [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
    "100bp DNA Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}

def clean_sequence(seq):
    if not seq: return ""
    return re.sub(r'[^a-zA-Z]', '', seq).upper()

def processar_upload(input_data):
    try:
        # Suporte a upload direto (BytesIO) ou arquivo do Streamlit
        nome_arquivo = input_data.name
        nome_sugerido = nome_arquivo.rsplit('.', 1)[0]
        
        if nome_arquivo.lower().endswith('.dna'):
            try:
                bytes_io = BytesIO(input_data.getvalue())
                record = SeqIO.read(bytes_io, "snapgene")
                return nome_sugerido, str(record.seq).upper()
            except Exception as e: return "Erro", f"Erro .dna: {str(e)}"
        
        bytes_data = input_data.getvalue()
        try: conteudo = bytes_data.decode("utf-8")
        except: conteudo = bytes_data.decode("latin-1")
        
        if ">" in conteudo:
            try:
                iterator = SeqIO.parse(StringIO(conteudo), "fasta")
                record = next(iterator)
                return record.id, str(record.seq).upper()
            except: pass
            
        seq_limpa = "".join([l.strip() for l in conteudo.splitlines() if not l.startswith(">")])
        seq_final = clean_sequence(seq_limpa)
        if len(seq_final) > 0 and any(c not in "ATGCNRYKMSWBDHV" for c in seq_final[:100]):
             return "Erro", "Arquivo inválido."
        return nome_sugerido, seq_final
    except Exception as e: return "Erro", str(e)

def processar_texto_manual(texto):
    if ">" in texto:
        try:
            iterator = SeqIO.parse(StringIO(texto), "fasta")
            record = next(iterator)
            return record.id, str(record.seq).upper()
        except: pass
    return "Seq Manual", clean_sequence(texto)

def calcular_digestao(sequencia, enzimas, eh_circular):
    if not sequencia or sequencia.startswith("Erro"): return []
    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)
    if eh_circular and not enzimas: return [(tamanho_total * 1.4, "Nicked", tamanho_total), (tamanho_total * 0.7, "Supercoiled", tamanho_total)]
    if not enzimas: return [(tamanho_total, "Linear", tamanho_total)]
    
    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))
    if not locais: return [(tamanho_total, "Uncut", tamanho_total)]
        
    fragmentos = []
    if not eh_circular:
        prev = 0
        for cut in locais:
            fragmentos.append(cut - prev)
            prev = cut
        fragmentos.append(tamanho_total - prev)
    else:
        if len(locais) == 1: fragmentos.append(tamanho_total)
        else:
            for i in range(len(locais)-1): fragmentos.append(locais[i+1] - locais[i])
            fragmentos.append((tamanho_total - locais[-1]) + locais[0])
    return [(frag, "Fragmento", frag) for frag in sorted(fragmentos, reverse=True)]

def smart_pcr_search(template, fwd_raw, rev_raw, eh_circular):
    fwd = clean_sequence(fwd_raw)
    rev = clean_sequence(rev_raw)
    template = template.upper()
    diag = {'fwd_found': False, 'rev_found': False, 'products': 0}
    if len(fwd) < 10 or len(rev) < 10: return [], diag

    SEED = 12
    fwd_seed = fwd[-SEED:]
    fwd_matches = [m.start() for m in re.finditer(fwd_seed, template)]
    if fwd_matches: diag['fwd_found'] = True
    
    rev_matches = []
    rev_seed_std = rev[-SEED:]
    rev_rc_std = str(Seq(rev_seed_std).reverse_complement())
    matches_std = [m.start() for m in re.finditer(rev_rc_std, template)]
    
    rev_seed_inv = rev[:SEED]
    rev_compl_inv = str(Seq(rev_seed_inv).complement())
    matches_inv = [m.start() for m in re.finditer(rev_compl_inv, template)]
    
    if len(matches_inv) > 0: rev_matches = matches_inv
    elif len(matches_std) > 0: rev_matches = matches_std
    if rev_matches: diag['rev_found'] = True
    
    produtos = []
    for f_pos in fwd_matches:
        f_end = f_pos + len(fwd_seed)
        for r_pos in rev_matches:
            if r_pos > f_pos:
                dist = r_pos - f_end
                if dist >= 0: produtos.append(len(fwd) + len(rev) + dist)
            elif eh_circular and r_pos < f_pos:
                dist = (len(template) - f_end) + r_pos
                produtos.append(len(fwd) + len(rev) + dist)
    diag['products'] = len(produtos)
    return [(p, "PCR", p) for p in sorted(produtos, reverse=True)], diag

# --- 6. INTERFACE ---

if 'lang' not in st.session_state: st.session_state.lang = "PT"

with st.sidebar:
    st.markdown("""<div style="text-align: left; margin-bottom: 10px;"><h1 style="color: #0F766E; margin:0; font-size:24px;">BioSpark</h1><p style="font-size: 10px; color: #0F766E; margin-top:-2px;">STUDIO v2.0</p></div>""", unsafe_allow_html=True)
    st.markdown("---")
    lang = st.session_state.lang
    
    # Abas da Sidebar: Configuração e Caderno
    tab_conf, tab_note = st.tabs(["⚙️ " + TEXTS["sidebar_config"][lang], TEXTS["notebook_title"][lang]])
    
    with tab_conf:
        num_pocos = st.slider(TEXTS["sidebar_wells"][lang], 1, 15, 4)
        agarose = st.slider(TEXTS["sidebar_agarose"][lang], 0.5, 2.0, 1.0, 0.1)
        st.caption(TEXTS["sidebar_visual"][lang])
        estilo_gel = st.selectbox(TEXTS["sidebar_theme"][lang], ["Profissional (Dark P&B)", "Publicação (Light P&B)", "Neon (Verde/Laranja)"])
        
        st.markdown("---")
        with st.expander(f"ℹ️ {TEXTS['guide_title'][lang]}"): st.markdown(TEXTS["guide_content"][lang])
        new_lang = st.selectbox("Lang", ["Português", "English"], index=0 if lang=="PT" else 1, label_visibility="collapsed")
        if (new_lang == "Português" and lang != "PT") or (new_lang == "English" and lang != "EN"):
            st.session_state.lang = "PT" if new_lang == "Português" else "EN"
            st.rerun()

    with tab_note:
        st.session_state.lab_notebook = st.text_area("Notas:", st.session_state.lab_notebook, height=300, placeholder="Anote seus primers, conclusões e observações aqui...")

    # RODAPÉ LATERAL
    st.markdown(f"""
    <div class="sidebar-footer">
        {TEXTS['created_by'][lang]} <strong>Elton Ostetti</strong>
        <a class="bug-report" href="mailto:e.ostetti.proppg@proppg.butantan.gov.br?subject=Bug%20Report%20BioSpark">{TEXTS['report_bug'][lang]}</a>
        <br>
        <span style="opacity:0.9;">FAPESP • USP • Instituto Butantan</span>
    </div>
    """, unsafe_allow_html=True)

# --- MAIN ---
st.markdown(f"# {TEXTS['header_title'][lang]}")
st.markdown(TEXTS["header_sub"][lang])

# --- ÁREA DE UPLOAD EM LOTE INTELIGENTE ---
with st.expander(TEXTS['batch_label'][lang], expanded=False):
    # O widget key='batch_uploader' salva a lista de arquivos no session_state automaticamente
    st.file_uploader("Arraste arquivos (.dna, .fasta) para preencher os poços automaticamente", 
                     type=['dna', 'fasta', 'txt'], 
                     accept_multiple_files=True, 
                     key="batch_uploader", 
                     on_change=handle_batch_upload)

st.markdown(" ")

cols = st.columns(4)
relatorio = []
plot_data = []
labels_x = []

for i in range(num_pocos):
    with cols[i % 4]:
        with st.expander(f"🔹 {TEXTS['well_title'][lang]} {i+1}", expanded=(i==0)):
            tipo_disp = st.radio("Tipo", [TEXTS['opt_sample'][lang], TEXTS['opt_pcr'][lang], TEXTS['opt_ladder'][lang]], horizontal=True, label_visibility="collapsed", key=f"t{i}")
            
            if tipo_disp == TEXTS['opt_ladder'][lang]:
                lad = st.selectbox("Ladder", list(LADDERS.keys()), key=f"l{i}")
                plot_data.append([(t, "Ladder", t) for t in LADDERS[lad]])
                lbl = st.text_input(TEXTS['label_gel'][lang], "M", key=f"lb{i}")
                labels_x.append(lbl)
                relatorio.append({"Poço": i+1, "Tipo": "Ladder", "Detalhes": lad, "Bandas": str(LADDERS[lad])})
            
            else:
                tab1, tab2 = st.tabs([TEXTS['tab_file'][lang], TEXTS['tab_text'][lang]])
                seq, name = "", ""
                
                # UPLOAD INDIVIDUAL
                with tab1:
                    f = st.file_uploader(TEXTS['upload_label'][lang], type=['dna','fasta','txt'], key=f"u{i}", label_visibility="collapsed")
                    if f: name, seq = processar_upload(f)
                
                # EDITOR DE TEXTO (Conectado ao Session State para receber o Batch Upload)
                with tab2:
                    # key=f"tx_{i}" conecta este widget ao valor definido em handle_batch_upload
                    txt = st.text_area(TEXTS['paste_label'][lang], height=70, key=f"tx_{i}", label_visibility="collapsed", placeholder="ATGC...")
                    if txt and not seq: name, seq = processar_texto_manual(txt)
                
                st.markdown("---")
                # RÓTULO (Conectado ao Session State)
                val_rotulo = st.session_state.get(f"lbl_{i}", name[:15] if name else str(i+1))
                label = st.text_input(TEXTS['label_gel'][lang], val_rotulo, key=f"lbl_{i}")
                labels_x.append(label)
                
                if tipo_disp == TEXTS['opt_sample'][lang]: # Digestão
                    c = st.checkbox(TEXTS['check_circular'][lang], True, key=f"c{i}")
                    e = st.multiselect(TEXTS['sel_enzymes'][lang], TODAS_ENZIMAS, key=f"e{i}")
                    if seq:
                        res = calcular_digestao(seq, e, c)
                        plot_data.append(res)
                        relatorio.append({"Poço": i+1, "Tipo": "Digestão", "Bandas": ";".join([str(int(x[0])) for x in res])})
                    else: 
                        plot_data.append([])
                        relatorio.append({"Poço": i+1, "Tipo": "Vazio"})
                        
                else: # PCR
                    fwd = st.text_input(TEXTS['pcr_fwd'][lang], key=f"fw{i}")
                    rev = st.text_input(TEXTS['pcr_rev'][lang], key=f"rv{i}")
                    c = st.checkbox(TEXTS['check_circular'][lang], False, key=f"cp{i}")
                    
                    if seq and fwd and rev:
                        res, diag = smart_pcr_search(seq, fwd, rev, c)
                        plot_data.append(res)
                        if not diag['fwd_found']: st.markdown(f"<p class='error-text'>{TEXTS['diag_fwd_fail'][lang]}</p>", unsafe_allow_html=True)
                        if not diag['rev_found']: st.markdown(f"<p class='error-text'>{TEXTS['diag_rev_fail'][lang]}</p>", unsafe_allow_html=True)
                        if diag['products'] > 1: st.markdown(f"<p class='warning-text'>{TEXTS['warn_multiple'][lang]}</p>", unsafe_allow_html=True)
                        elif diag['products'] == 0 and diag['fwd_found'] and diag['rev_found']: st.warning(TEXTS['warn_no_product'][lang])
                        relatorio.append({"Poço": i+1, "Tipo": "PCR", "Bandas": ";".join([str(int(x[0])) for x in res])})
                    else:
                        plot_data.append([])
                        relatorio.append({"Poço": i+1, "Tipo": "Vazio"})

# --- PLOTAGEM ---
st.markdown(" ")
if any(plot_data):
    if "Neon" in estilo_gel: bg, txt, c_samp, c_lad = '#111827', 'white', '#00ff41', '#ff9900'
    elif "Profissional" in estilo_gel: bg, txt, c_samp, c_lad = '#000000', 'white', 'white', 'white'
    else: bg, txt, c_samp, c_lad = 'white', 'black', 'black', 'black'

    min_view = 20
    max_view = 25000 / (agarose * 0.8)
    max_range = max(num_pocos, 15) + 0.5

    fig = go.Figure()
    for i, bands in enumerate(plot_data):
        x = i + 1
        is_ladder = (relatorio[i].get("Tipo") == "Ladder")
        color = c_lad if is_ladder else c_samp
        
        for (size, type, real) in bands:
            if size < (min_view * 0.9) or size > (max_view * 1.1): continue
            
            width = 2; opacity = 0.8
            if is_ladder:
                if size in [500, 1000, 3000]: width=5; opacity=1.0
            else:
                if type == "Supercoiled": width=4; opacity=0.7
                elif type == "PCR": width=3; opacity=0.9
            
            # ESTILO PALITO COM MARCADOR (RESTAURADO DA V1.0)
            fig.add_trace(go.Scatter(
                x=[x-0.28, x+0.28], y=[size, size], mode='lines+markers',
                line=dict(color=color, width=width), 
                marker=dict(color=color, size=width, symbol='circle'),
                opacity=opacity,
                hoverinfo='text', hovertext=f"{int(size)} pb",
                showlegend=False
            ))
            
            if is_ladder:
                fig.add_trace(go.Scatter(x=[x-0.4], y=[size], mode="text", text=[str(size)], textfont=dict(color=txt, size=9), showlegend=False))

    fig.update_layout(
        plot_bgcolor=bg, paper_bgcolor=bg, height=600,
        margin=dict(t=30, b=30, l=50, r=30),
        xaxis=dict(tickvals=list(range(1, num_pocos+1)), ticktext=labels_x, showgrid=False, tickfont=dict(color=txt), range=[0.5, max_range]),
        yaxis=dict(type='log', range=[math.log10(min_view), math.log10(max_view)], showgrid=False, showticklabels=False)
    )
    st.plotly_chart(fig, use_container_width=True)
    
    df = pd.DataFrame(relatorio)
    csv = df.to_csv(index=False).encode('utf-8')
    with st.expander(f"📥 {TEXTS['export_expander'][lang]}"):
        st.download_button(TEXTS['btn_download'][lang], csv, "gel_results.csv", "text/csv")
else:
    st.info(TEXTS['empty_msg'][lang])

st.markdown("""<div class="footer"><p><b>BioSpark</b></p></div>""", unsafe_allow_html=True)
