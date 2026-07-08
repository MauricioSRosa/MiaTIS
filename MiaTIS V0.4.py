import os
import re
import json
import pandas as pd
import requests
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from concurrent.futures import ThreadPoolExecutor
from collections import Counter

# ==========================================
# 1. PARSERS AVANÇADOS (HIERARQUIA BRITE E PATHWAY)
# ==========================================

KEGG_HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
}

def buscar_reacao_kegg(reaction_id):
    try:
        url = f"https://rest.kegg.jp/get/{reaction_id}"
        resp = requests.get(url, headers=KEGG_HEADERS, timeout=10)
        if resp.status_code != 200: return None
        match = re.search(r"^EQUATION\s+(.+)$", resp.text, re.MULTILINE)
        if not match: return None
        partes = re.split(r'<=>', match.group(1))
        if len(partes) != 2: return None
        return {
            "substrates": re.findall(r'C\d{5}', partes[0]),
            "products": re.findall(r'C\d{5}', partes[1])
        }
    except: return None

def buscar_compostos_kegg_lote(lista_ids):
    if not lista_ids: return {}
    resultados_bloco = {}
    tamanho_lote = 20
    
    for i in range(0, len(lista_ids), tamanho_lote):
        sub_lista = lista_ids[i:i+tamanho_lote]
        url = f"https://rest.kegg.jp/get/{'+'.join(sub_lista)}"
        try:
            resp = requests.get(url, headers=KEGG_HEADERS, timeout=15)
            if resp.status_code != 200: continue
            entradas = resp.text.split('///')
            for entrada in entradas:
                if not entrada.strip(): continue
                linhas = entrada.split('\n')
                comp_id = None
                nomes = []
                pathways = []
                brite_lines = []
                current_key = None
                
                for linha in linhas:
                    if not linha.strip(): continue
                    if not linha.startswith(' '):
                        current_key = linha[:12].strip()
                        conteudo = linha[12:]
                    else:
                        conteudo = linha[12:]
                        
                    if current_key == 'ENTRY':
                        partes = conteudo.split()
                        if partes: comp_id = partes[0].strip()
                    elif current_key == 'NAME':
                        nomes.append(conteudo.strip())
                    elif current_key == 'PATHWAY':
                        pathways.append(conteudo.strip())
                    elif current_key == 'BRITE':
                        brite_lines.append(conteudo)
                        
                if comp_id:
                    nome_completo = " ".join(nomes)
                    lista_nomes = [n.strip() for n in nome_completo.split(';') if n.strip()]
                    
                    # --- RECONSTRUTOR DE HIERARQUIA BRITE BASEADO EM INDENTAÇÃO ---
                    brite_paths = []
                    current_tree = {}
                    for b_line in brite_lines:
                        # Conta os espaços adicionais de recuo a partir da coluna de conteúdo (index 12+)
                        leading_spaces = len(b_line) - len(b_line.lstrip(' '))
                        txt_limpo = b_line.strip()
                        if not txt_limpo: continue
                        
                        level = leading_spaces
                        current_tree[level] = txt_limpo
                        
                        # Limpa ramos residuais de ramificações paralelas mais profundas
                        for k in list(current_tree.keys()):
                            if k > level: del current_tree[k]
                            
                        # Identifica o nó folha do composto (usa os dígitos numéricos para bater Cxxxxx ou Dxxxxx)
                        if len(comp_id) >= 6 and comp_id[1:] in txt_limpo:
                            path = [current_tree[l] for l in sorted(current_tree.keys()) if l < level]
                            if path: brite_paths.append(path)
                            
                    resultados_bloco[comp_id] = {
                        "names": " / ".join(lista_nomes) if lista_nomes else comp_id,
                        "pathways": pathways,
                        "brite_hierarchies": brite_paths
                    }
        except: continue
    return resultados_bloco

# ==========================================
# 2. INTERFACE E GERENCIADOR DE DIRETÓRIOS
# ==========================================

class TnSeqMetabolicApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Tn-Seq Metabolic Pipeline - Análise Hierárquica Estratificada V5")
        self.geometry("950x800")
        
        self.master_cache_dir = os.path.join(os.getcwd(), "Banco_Global_KEGG")
        os.makedirs(self.master_cache_dir, exist_ok=True)
        self.rxn_cache_file = os.path.join(self.master_cache_dir, "master_reactions.json")
        self.comp_cache_file = os.path.join(self.master_cache_dir, "master_compounds.json")
        
        self.var_projeto_dir = tk.StringVar()
        self.var_transit_teste = tk.StringVar()
        self.var_transit_controle = tk.StringVar()
        self.var_eggnog_bac = tk.StringVar()
        self.var_eggnog_hospedeiro = tk.StringVar()
        self.var_corte_hubs = tk.IntVar(value=50)
        
        self.montar_interface()

    def montar_interface(self):
        abas = ttk.Notebook(self)
        abas.pack(expand=True, fill='both', padx=10, pady=5)
        
        aba1 = ttk.Frame(abas)
        abas.add(aba1, text="1. Parâmetros de Entrada")
        
        f_proj = ttk.LabelFrame(aba1, text="Diretório de Destino (Resultados)")
        f_proj.pack(fill='x', padx=10, pady=5)
        ttk.Entry(f_proj, textvariable=self.var_projeto_dir, width=75).pack(side='left', padx=5, pady=5)
        ttk.Button(f_proj, text="Procurar...", command=lambda: self.abrir_diretorio(self.var_projeto_dir)).pack(side='left', padx=5)
        
        f_trans = ttk.LabelFrame(aba1, text="Arquivos TRANSIT (.genes.txt)")
        f_trans.pack(fill='x', padx=10, pady=5)
        ttk.Label(f_trans, text="Interação (In Vivo):").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_trans, textvariable=self.var_transit_teste, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(f_trans, text="Abrir", command=lambda: self.abrir_arquivo(self.var_transit_teste)).grid(row=0, column=2, padx=5)
        ttk.Label(f_trans, text="Controle (In Vitro):").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_trans, textvariable=self.var_transit_controle, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(f_trans, text="Abrir", command=lambda: self.abrir_arquivo(self.var_transit_controle)).grid(row=1, column=2, padx=5)
        
        f_egg = ttk.LabelFrame(aba1, text="Anotações Funcionais EggNOG (.tabular)")
        f_egg.pack(fill='x', padx=10, pady=5)
        ttk.Label(f_egg, text="Bactéria:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_egg, textvariable=self.var_eggnog_bac, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(f_egg, text="Abrir", command=lambda: self.abrir_arquivo(self.var_eggnog_bac)).grid(row=0, column=2, padx=5)
        ttk.Label(f_egg, text="Hospedeiro (Planta):").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_egg, textvariable=self.var_eggnog_hospedeiro, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(f_egg, text="Abrir", command=lambda: self.abrir_arquivo(self.var_eggnog_hospedeiro)).grid(row=1, column=2, padx=5)
        
        aba2 = ttk.Frame(abas)
        abas.add(aba2, text="2. Processamento")
        f_param = ttk.LabelFrame(aba2, text="Filtros Topológicos")
        f_param.pack(fill='x', padx=10, pady=5)
        ttk.Label(f_param, text="Limite de Conexões para Hubs:").pack(side='left', padx=5, pady=5)
        ttk.Entry(f_param, textvariable=self.var_corte_hubs, width=8).pack(side='left', padx=5, pady=5)
        ttk.Button(aba2, text="RODAR MAPEAMENTO MULTI-ESTRATIFICADO V5", command=self.rodar_pipeline).pack(pady=15, ipady=5)
        
        f_log = ttk.LabelFrame(self, text="Console de Acompanhamento")
        f_log.pack(fill='both', expand=True, padx=10, pady=5)
        self.txt_log = tk.Text(f_log, height=14, bg='black', fg='#00FF00', font=('Consolas', 10))
        self.txt_log.pack(fill='both', expand=True, padx=5, pady=5)

    def log(self, msg):
        self.txt_log.insert(tk.END, msg + "\n")
        self.txt_log.see(tk.END)
        self.update_idletasks()

    def abrir_arquivo(self, var):
        p = filedialog.askopenfilename()
        if p: var.set(p)

    def abrir_diretorio(self, var):
        p = filedialog.askdirectory()
        if p: var.set(p)

    def carregar_banco_master(self):
        rxn = json.load(open(self.rxn_cache_file, 'r')) if os.path.exists(self.rxn_cache_file) else {}
        comp = json.load(open(self.comp_cache_file, 'r')) if os.path.exists(self.comp_cache_file) else {}
        return rxn, comp

    def salvar_banco_master(self, rxn, comp):
        with open(self.rxn_cache_file, 'w') as f: json.dump(rxn, f)
        with open(self.comp_cache_file, 'w') as f: json.dump(comp, f)

    def carregar_dados_essencialidade(self, out_dir):
        self.log(" -> Carregando arquivo TRANSIT (.genes.txt)...")
        df_teste = pd.read_csv(self.var_transit_teste.get(), sep="\t", skiprows=4)
        
        col_gene = 'gene' if 'gene' in df_teste.columns else ('#ORF' if '#ORF' in df_teste.columns else df_teste.columns[0])
        df_teste.rename(columns={col_gene: 'ID_GENE'}, inplace=True)
        
        # TRATAMENTO DE ERRO: Garante que valores vazios ou NaN na coluna call virem '--'
        if 'call' in df_teste.columns:
            df_teste['call'] = df_teste['call'].fillna("--")
        
        if self.var_transit_controle.get():
            df_ctrl = pd.read_csv(self.var_transit_controle.get(), sep="\t", skiprows=4)
            col_ctrl = 'gene' if 'gene' in df_ctrl.columns else ('#ORF' if '#ORF' in df_ctrl.columns else df_ctrl.columns[0])
            df_ctrl.rename(columns={col_ctrl: 'ID_GENE'}, inplace=True)
            
            if 'call' in df_ctrl.columns:
                df_ctrl['call'] = df_ctrl['call'].fillna("--")
            
            df_m = pd.merge(df_ctrl, df_teste, on="ID_GENE", suffixes=("_controle", "_condicao"))
            df_m = df_m[~((df_m["call_controle"] == "ES") & (df_m["call_condicao"] == "ES"))]
            df_m.to_csv(os.path.join(out_dir, "Essencialidade_Filtro_Controle.txt"), sep="\t", index=False)
            return dict(zip(df_m['ID_GENE'], df_m['call_condicao'].fillna("--")))
        else:
            # Garante que o dicionário final não leve valores float 'nan'
            return dict(zip(df_teste['ID_GENE'], df_teste['call'].fillna("--")))

    def extrair_rede_para_csv(self, arquivo_eggnog, cache_rxn, cache_comp, caminho_csv, organism_label):
        self.log(f" -> Processando cache estruturado para {organism_label}...")
        with open(arquivo_eggnog, 'r', encoding='utf-8') as f:
            lin = f.readline()
            cols = [c.strip() for c in lin.lstrip('#').strip().split('\t')] if lin.startswith('#') else None
        
        df = pd.read_csv(arquivo_eggnog, sep="\t", comment="#", names=cols, skiprows=1) if cols else pd.read_csv(arquivo_eggnog, sep="\t")
        col_query = 'query' if 'query' in df.columns else df.columns[0]
        df_limpo = df[[col_query, 'KEGG_Reaction']].dropna()
        
        reacoes_alvo = set()
        for rxns in df_limpo['KEGG_Reaction'].astype(str):
            for r in rxns.split(','):
                if re.match(r'^R\d{5}$', r.strip()): reacoes_alvo.add(r.strip())
                
        novas_rxn = [r for r in reacoes_alvo if r not in cache_rxn]
        if novas_rxn:
            self.log(f"   + Buscando {len(novas_rxn)} equações estruturais no KEGG...")
            with ThreadPoolExecutor(max_workers=4) as ex:
                res = list(ex.map(buscar_reacao_kegg, novas_rxn))
                for r_id, r_dados in zip(novas_rxn, res):
                    cache_rxn[r_id] = r_dados if r_dados else {"substrates": [], "products": []}

        linhas_rede = []
        compostos_detectados = set()
        
        for _, row in df_limpo.iterrows():
            enzima = str(row.iloc[0]).strip()
            reacoes = [r.strip() for r in str(row['KEGG_Reaction']).split(',')]
            for r in reacoes:
                if r in cache_rxn:
                    subs = cache_rxn[r].get('substrates', [])
                    prods = cache_rxn[r].get('products', [])
                    if subs or prods:
                        compostos_detectados.update(subs + prods)
                        linhas_rede.append({
                            "Substratos": ";".join(subs) if subs else "Nenhum",
                            "Enzima": enzima,
                            "Produtos": ";".join(prods) if prods else "Nenhum"
                        })
                        
        df_rede = pd.DataFrame(linhas_rede)
        df_rede.to_csv(caminho_csv, index=False)
        
        # Garante o re-download completo caso a estrutura salva no JSON local seja de uma versão anterior sem hierarquias
        novos_comp = [c for c in compostos_detectados if c not in cache_comp or "brite_hierarchies" not in cache_comp[c] or "pathways" not in cache_comp[c]]
        if novos_comp:
            self.log(f"   + Atualizando metadados estruturais de {len(novos_comp)} compostos via KEGG API...")
            res_lotes = buscar_compostos_kegg_lote(novos_comp)
            for c_id, d in res_lotes.items(): cache_comp[c_id] = d
            for c_id in novos_comp:
                if c_id not in cache_comp: 
                    cache_comp[c_id] = {"names": c_id, "pathways": [], "brite_hierarchies": []}
                
        return df_rede, cache_rxn, cache_comp

    # ==========================================
    # 3. CONSTRUÇÃO DA MATRIZ MESTRE DE CLASSIFICAÇÃO
    # ==========================================

    def executar_tripla_classificacao(self, df_bac, mapa_essencialidade, cache_comp, out_dir):
        self.log("\n[CONSTRUINDO ARQUITETURA DE DECISÃO DOS METABÓLITOS]")
        
        todos_mets_brutos = []
        for _, r in df_bac.iterrows():
            if r['Substratos'] != "Nenhum": todos_mets_brutos.extend(r['Substratos'].split(';'))
            if r['Produtos'] != "Nenhum": todos_mets_brutos.extend(r['Produtos'].split(';'))
        
        contagem_hubs = Counter(todos_mets_brutos)
        limite_corte = self.var_corte_hubs.get()
        hubs_removidos = {m for m, c in contagem_hubs.items() if c > limite_corte}
        self.log(f" -> Filtro de Hubs: {len(hubs_removidos)} metabólitos ubíquos desconsiderados (> {limite_corte} conexões).")

        enzimas_da_rede = set(df_bac['Enzima'].unique())
        linhas_genes_csv = []
        for enz in enzimas_da_rede:
            classe_ess = mapa_essencialidade.get(enz, "--")
            linhas_genes_csv.append({"Gene": enz, "Essencialidade": classe_ess})
        df_genes_mapeados = pd.DataFrame(linhas_genes_csv)
        df_genes_mapeados.to_csv(os.path.join(out_dir, "Mapeamento_Essencialidade_Genes.csv"), index=False)
        mapa_genes_final = dict(zip(df_genes_mapeados['Gene'], df_genes_mapeados['Essencialidade']))

        set_substratos = set()
        set_produtos = set()
        for _, row in df_bac.iterrows():
            if row['Substratos'] != "Nenhum":
                set_substratos.update([s for s in row['Substratos'].split(';') if s not in hubs_removidos])
            if row['Produtos'] != "Nenhum":
                set_produtos.update([p for p in row['Produtos'].split(';') if p not in hubs_removidos])
                
        todos_metabolitos = set_substratos.union(set_produtos)
        hierarquia = {"ES": 4, "GD": 3, "NE": 2, "GA": 1, "--": 0, "SS": -1, "PP": -1}
        inv_hierarquia = {4: "ES", 3: "GD", 2: "NE", 1: "GA", 0: "--"}

        linhas_metabolitos_master = []
        for met in todos_metabolitos:
            if met in set_substratos and met not in set_produtos:
                class_topologica = "S"
            elif met in set_produtos and met not in set_substratos:
                class_topologica = "P"
            else:
                class_topologica = "M"

            if class_topologica == "S":
                class_produtores = "SS"
            else:
                enzimas_produtoras = df_bac[df_bac['Produtos'].str.contains(met, na=False)]['Enzima'].unique()
                valores_ess = []
                for e in enzimas_produtoras:
                    val = mapa_genes_final.get(e, "--")
                    if pd.isna(val) or val not in hierarquia:
                        val = "--"
                    valores_ess.append(hierarquia[val])
                class_produtores = inv_hierarquia[max(valores_ess)] if valores_ess else "--"

            if class_topologica == "P":
                class_consumidores = "PP"
            else:
                enzimas_consumidoras = df_bac[df_bac['Substratos'].str.contains(met, na=False)]['Enzima'].unique()
                valores_ess = []
                for e in enzimas_consumidoras:
                    val = mapa_genes_final.get(e, "--")
                    if pd.isna(val) or val not in hierarquia:
                        val = "--"
                    valores_ess.append(hierarquia[val])
                class_consumidores = inv_hierarquia[max(valores_ess)] if valores_ess else "--"

            info_comp = cache_comp.get(met, {})
            linhas_metabolitos_master.append({
                "Metabolito": met,
                "Nome": info_comp.get('names', met),
                "Classificacao_Topologica": class_topologica,
                "Essencialidade_Produtores": class_produtores,
                "Essencialidade_Consumidores": class_consumidores
            })

        df_master_metabolitos = pd.DataFrame(linhas_metabolitos_master)
        df_master_metabolitos.to_csv(os.path.join(out_dir, "Classificacao_Metabolitos_Master.csv"), index=False)
        return df_master_metabolitos, hubs_removidos, mapa_genes_final

    # ==========================================
    # 4. EXPORTAÇÃO E ENRIQUECIMENTO DOS 6 SUBGRUPOS
    # ==========================================

    def exportar_enriquecimento_grupo(self, lista_compounds, cache_comp, prefixo, out_dir):
        """Gera os relatórios de nomes, hierarquias BRITE (por nível) e PATHWAY para um subgrupo."""
        pasta_grupo = os.path.join(out_dir, "Estratificacao_Subgrupos")
        os.makedirs(pasta_grupo, exist_ok=True)
        
        # 1. Exportação da Lista de Nomes
        linhas_comp = []
        for c_id in lista_compounds:
            info = cache_comp.get(c_id, {})
            linhas_comp.append({"Metabolito": c_id, "Nome": info.get("names", c_id)})
        df_lista = pd.DataFrame(linhas_comp)
        df_lista.to_csv(os.path.join(pasta_grupo, f"Lista_Compostos_{prefixo}.csv"), index=False, encoding='utf-8')
        
        if not lista_compounds: return
        
        # 2. Processamento Hierárquico do BRITE
        level_counts = {}
        pathway_counts = Counter()
        
        for c_id in lista_compounds:
            info = cache_comp.get(c_id, {})
            for pw in info.get("pathways", []):
                pathway_counts[pw] += 1
                
            for path in info.get("brite_hierarchies", []):
                for level_idx, categoria in enumerate(path):
                    if level_idx not in level_counts:
                        level_counts[level_idx] = Counter()
                    level_counts[level_idx][categoria] += 1
                    
        # Salva Enriquecimento BRITE por Nível
        linhas_brite = []
        for lvl in sorted(level_counts.keys()):
            for cat, count in level_counts[lvl].most_common():
                linhas_brite.append({
                    "Nivel_Hierarquia": f"Nivel_{lvl}",
                    "Categoria_BRITE": cat,
                    "Ocorrencias": count
                })
        df_brite = pd.DataFrame(linhas_brite)
        df_brite.to_csv(os.path.join(pasta_grupo, f"Enriquecimento_BRITE_{prefixo}.csv"), index=False, encoding='utf-8')
        
        # Salva Enriquecimento PATHWAY (Plano)
        linhas_pw = []
        for pw, count in pathway_counts.most_common():
            linhas_pw.append({"Via_Metabolica": pw, "Ocorrencias": count})
        df_pw = pd.DataFrame(linhas_pw)
        df_pw.to_csv(os.path.join(pasta_grupo, f"Enriquecimento_PATHWAY_{prefixo}.csv"), index=False, encoding='utf-8')

    def processar_gargalos_finais(self, df_master_met, df_bac, hubs_removidos, mapa_genes, produtos_hospedeiro, cache_comp, out_dir):
        self.log(" -> Isolando e classificando os pontos de gargalo biológico...")
        
        linhas_gargalos = []
        df_filtrado_gargalo = df_master_met[
            (df_master_met['Essencialidade_Consumidores'].isin(['ES', 'GD'])) &
            (df_master_met['Essencialidade_Produtores'].isin(['NE', 'GA', 'SS', '--']))
        ]

        for _, row in df_filtrado_gargalo.iterrows():
            met = row['Metabolito']
            class_prod = row['Essencialidade_Produtores']
            
            enz_consumidoras = df_bac[df_bac['Substratos'].str.contains(met, na=False)]['Enzima'].unique()
            enz_vulneraveis = [e for e in enz_consumidoras if mapa_genes.get(e, "--") in ['ES', 'GD']]
            
            if class_prod == "SS":
                tipo_gargalo = "Base Externa (Substrato não sintetizado)"
            else:
                enz_produtoras = df_bac[df_bac['Produtos'].str.contains(met, na=False)]['Enzima'].unique()
                precursores_vulneraveis = False
                for ep in enz_produtoras:
                    subs_da_ep = df_bac[df_bac['Enzima'] == ep]['Substratos'].values
                    for sub_lista in subs_da_ep:
                        if sub_lista == "Nenhum": continue
                        for s in sub_lista.split(';'):
                            if s in hubs_removidos: continue
                            if df_master_met[df_master_met['Metabolito'] == s]['Essencialidade_Produtores'].isin(['ES', 'GD']).any():
                                precursores_vulneraveis = True
                                
                tipo_gargalo = "Espremida (Importante -> Fraco -> Importante)" if precursores_vulneraveis else "Base Interna (Fraco -> Fraco -> Importante)"

            esta_na_planta = met in produtos_hospedeiro
            linhas_gargalos.append({
                "Metabolito_Interface": met,
                "Nome_Metabolito": cache_comp.get(met, {}).get('names', met),
                "Classificacao_Topologica": row['Classificacao_Topologica'],
                "Tipo_Gargalo": tipo_gargalo,
                "Presente_No_Hospedeiro": "SIM" if esta_na_planta else "NÃO",
                "Enzimas_Importantes_Afetadas": ";".join(enz_vulneraveis)
            })

        df_gargalos_total = pd.DataFrame(linhas_gargalos)
        df_gargalos_total.to_csv(os.path.join(out_dir, "Mapeamento_Gargalos_Metabolicos.csv"), index=False, encoding='utf-8')

        # --- SEPARAÇÃO DAS ANÁLISES NOS 6 SUBGRUPOS REQUISITADOS ---
        gargalos_validos = df_gargalos_total[df_gargalos_total['Classificacao_Topologica'].isin(['S', 'M'])]
        
        subgrupos_matriz = {
            "Gargalos_Geral_SM_Presentes_Hospedeiro": gargalos_validos[gargalos_validos['Presente_No_Hospedeiro'] == "SIM"],
            "Gargalos_Geral_SM_Ausentes_Hospedeiro": gargalos_validos[gargalos_validos['Presente_No_Hospedeiro'] == "NÃO"],
            
            "Somente_S_Presentes_Hospedeiro": gargalos_validos[(gargalos_validos['Classificacao_Topologica'] == "S") & (gargalos_validos['Presente_No_Hospedeiro'] == "SIM")],
            "Somente_S_Ausentes_Hospedeiro": gargalos_validos[(gargalos_validos['Classificacao_Topologica'] == "S") & (gargalos_validos['Presente_No_Hospedeiro'] == "NÃO")],
            
            "Somente_M_Presentes_Hospedeiro": gargalos_validos[(gargalos_validos['Classificacao_Topologica'] == "M") & (gargalos_validos['Presente_No_Hospedeiro'] == "SIM")],
            "Somente_M_Ausentes_Hospedeiro": gargalos_validos[(gargalos_validos['Classificacao_Topologica'] == "M") & (gargalos_validos['Presente_No_Hospedeiro'] == "NÃO")]
        }

        self.log("\n[EXECUTANDO DESMISTIFICAÇÃO E ENRIQUECIMENTO POR SUBGRUPO]")
        relatorio_linhas = [
            "=========================================================",
            "  RELATÓRIO DE ESTRATIFICAÇÃO METABÓLICA TRIPLA (V5)     ",
            "=========================================================",
            f"Total Geral de Gargalos Analisados (S+M): {len(gargalos_validos)}\n"
        ]

        for nome_grupo, df_sub in subgrupos_matriz.items():
            lista_ids = df_sub['Metabolito_Interface'].tolist()
            self.exportar_enriquecimento_grupo(lista_ids, cache_comp, nome_grupo, out_dir)
            msg = f"  -> {nome_grupo}: {len(lista_ids)} compostos mapeados."
            self.log(msg)
            relatorio_linhas.append(msg)

        with open(os.path.join(out_dir, "Sumario_Estatistico_Gargalos.txt"), "w", encoding='utf-8') as f:
            f.write("\n".join(relatorio_linhas))

    # ==========================================
    # 5. EXECUÇÃO CENTRAL DA PIPELINE
    # ==========================================

    def rodar_pipeline(self):
        if not self.var_projeto_dir.get() or not self.var_transit_teste.get() or not self.var_eggnog_bac.get():
            messagebox.showerror("Erro", "Preencha os inputs obrigatórios antes de rodar.")
            return
            
        out_dir = os.path.join(self.var_projeto_dir.get(), "Pipeline_Outputs_V5")
        os.makedirs(out_dir, exist_ok=True)
        
        self.txt_log.delete('1.0', tk.END)
        self.log("[INICIANDO ENGINE V5 - DECODIFICADOR HIERÁRQUICO BRITE & PATHWAY]")

        rxn_cache, comp_cache = self.carregar_banco_master()
        mapa_essencialidade = self.carregar_dados_essencialidade(out_dir)

        df_bac_rede, rxn_cache, comp_cache = self.extrair_rede_para_csv(
            self.var_eggnog_bac.get(), rxn_cache, comp_cache, 
            os.path.join(out_dir, "Rede_Cache_Bacteria.csv"), "Bactéria"
        )
        
        produtos_hospedeiro = set()
        if self.var_eggnog_hospedeiro.get():
            df_hosp_rede, rxn_cache, comp_cache = self.extrair_rede_para_csv(
                self.var_eggnog_hospedeiro.get(), rxn_cache, comp_cache,
                os.path.join(out_dir, "Rede_Cache_Hospedeiro.csv"), "Hospedeiro"
            )
            for _, r in df_hosp_rede.iterrows():
                if r['Produtos'] != "Nenhum": produtos_hospedeiro.update(r['Produtos'].split(';'))

        self.salvar_banco_master(rxn_cache, comp_cache)

        df_master_met, hubs_removidos, mapa_genes = self.executar_tripla_classificacao(
            df_bac_rede, mapa_essencialidade, comp_cache, out_dir
        )

        self.processar_gargalos_finais(
            df_master_met, df_bac_rede, hubs_removidos, 
            mapa_genes, produtos_hospedeiro, comp_cache, out_dir
        )
        
        self.log("\n=== CONCLUÍDO! Verifique a pasta 'Pipeline_Outputs_V5/Estratificacao_Subgrupos' ===")
        messagebox.showinfo("Sucesso", "Análise Concluída! Subgrupos exportados com sucesso.")

if __name__ == "__main__":
    app = TnSeqMetabolicApp()
    app.mainloop()