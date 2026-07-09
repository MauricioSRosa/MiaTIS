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
# 1. Retrieve KEEG data
# ==========================================

KEGG_HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
}

def BucaReacoesKEGG(IDReacao):
    try:
        url = f"https://rest.kegg.jp/get/{IDReacao}"
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

def BuscarCompostosKeggLote(IDsLista):
    if not IDsLista: return {}
    DicionarioBlocoDeResultados = {}
    TamanhoLoteDeBusca = 20
    
    for i in range(0, len(IDsLista), TamanhoLoteDeBusca):
        ListaB = IDsLista[i:i+TamanhoLoteDeBusca]
        url = f"https://rest.kegg.jp/get/{'+'.join(ListaB)}"
        try:
            resp = requests.get(url, headers=KEGG_HEADERS, timeout=15)
            if resp.status_code != 200: continue
            entradas = resp.text.split('///')
            for entrada in entradas:
                if not entrada.strip(): continue
                linhas = entrada.split('\n')
                IDComposto = None
                Nomes = []
                RotaMetabolica = []
                BRITEMetabolitos = []
                ChaveAtual = None
                
                for linha in linhas:
                    if not linha.strip(): continue
                    if not linha.startswith(' '):
                        ChaveAtual = linha[:12].strip()
                        conteudo = linha[12:]
                    else:
                        conteudo = linha[12:]
                        
                    if ChaveAtual == 'ENTRY':
                        partes = conteudo.split()
                        if partes: IDComposto = partes[0].strip()
                    elif ChaveAtual == 'NAME':
                        Nomes.append(conteudo.strip())
                    elif ChaveAtual == 'PATHWAY':
                        RotaMetabolica.append(conteudo.strip())
                    elif ChaveAtual == 'BRITE':
                        BRITEMetabolitos.append(conteudo)
                        
                if IDComposto:
                    nome_completo = " ".join(Nomes)
                    lista_Nomes = [n.strip() for n in nome_completo.split(';') if n.strip()]
                    
                    brite_paths = []
                    current_tree = {}
                    for b_line in BRITEMetabolitos:
                        leading_spaces = len(b_line) - len(b_line.lstrip(' '))
                        txt_limpo = b_line.strip()
                        if not txt_limpo: continue
                        
                        level = leading_spaces
                        current_tree[level] = txt_limpo
                        
                        for k in list(current_tree.keys()):
                            if k > level: del current_tree[k]
                            
                        if len(IDComposto) >= 6 and IDComposto[1:] in txt_limpo:
                            path = [current_tree[l] for l in sorted(current_tree.keys()) if l < level]
                            if path: brite_paths.append(path)
                            
                    DicionarioBlocoDeResultados[IDComposto] = {
                        "names": " / ".join(lista_Nomes) if lista_Nomes else IDComposto,
                        "RotaMetabolica": RotaMetabolica,
                        "HierarquiaBRITE": brite_paths
                    }
        except: continue
    return DicionarioBlocoDeResultados

# ==========================================
# 2. Folder management and interface
# ==========================================

class TnSeqMetabolicApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Metabolic Interaction Analysis using TIS (MiaTIS) - V0.5")
        self.geometry("950x800")
        
        self.master_cache_dir = os.path.join(os.getcwd(), "KEGG_General_Directory")
        os.makedirs(self.master_cache_dir, exist_ok=True)
        self.rxn_cache_file = os.path.join(self.master_cache_dir, "master_reactions.json")
        self.comp_cache_file = os.path.join(self.master_cache_dir, "master_compounds.json")
        
        self.var_projeto_dir = tk.StringVar()
        self.var_transit_teste = tk.StringVar()
        self.var_transit_controle = tk.StringVar()
        self.var_eggnog_bac = tk.StringVar()
        self.var_eggnog_hospedeiro = tk.StringVar()
        self.var_corte_hubs = tk.IntVar(value=50)
        
        self.MontarInterface()

    def MontarInterface(self):
        abas = ttk.Notebook(self)
        abas.pack(expand=True, fill='both', padx=10, pady=5)
        
        aba1 = ttk.Frame(abas)
        abas.add(aba1, text="1. Input Data")
        
        f_proj = ttk.LabelFrame(aba1, text="Destination Directory (Output)")
        f_proj.pack(fill='x', padx=10, pady=5)
        ttk.Entry(f_proj, textvariable=self.var_projeto_dir, width=75).pack(side='left', padx=5, pady=5)
        ttk.Button(f_proj, text="Search ...", command=lambda: self.AbrirDiretorio(self.var_projeto_dir)).pack(side='left', padx=5)
        
        f_trans = ttk.LabelFrame(aba1, text="TRANSIT files (.genes.txt)")
        f_trans.pack(fill='x', padx=10, pady=5)
        ttk.Label(f_trans, text="Interaction (In Vivo):").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_trans, textvariable=self.var_transit_teste, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(f_trans, text="Open", command=lambda: self.AbrirArquivo(self.var_transit_teste)).grid(row=0, column=2, padx=5)
        ttk.Label(f_trans, text="Control (In Vitro):").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_trans, textvariable=self.var_transit_controle, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(f_trans, text="Open", command=lambda: self.AbrirArquivo(self.var_transit_controle)).grid(row=1, column=2, padx=5)
        
        f_egg = ttk.LabelFrame(aba1, text="EggNOG Functional Annotations (.tabular)")
        f_egg.pack(fill='x', padx=10, pady=5)
        ttk.Label(f_egg, text="Bacteria:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_egg, textvariable=self.var_eggnog_bac, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(f_egg, text="Open", command=lambda: self.AbrirArquivo(self.var_eggnog_bac)).grid(row=0, column=2, padx=5)
        ttk.Label(f_egg, text="Host:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(f_egg, textvariable=self.var_eggnog_hospedeiro, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(f_egg, text="Open", command=lambda: self.AbrirArquivo(self.var_eggnog_hospedeiro)).grid(row=1, column=2, padx=5)
        
        aba2 = ttk.Frame(abas)
        abas.add(aba2, text="2. Processing")
        f_param = ttk.LabelFrame(aba2, text="Topological Filters")
        f_param.pack(fill='x', padx=10, pady=5)
        ttk.Label(f_param, text="Connection Limit for Hubs:").pack(side='left', padx=5, pady=5)
        ttk.Entry(f_param, textvariable=self.var_corte_hubs, width=8).pack(side='left', padx=5, pady=5)
        ttk.Button(aba2, text="RUN MAPPING", command=self.rodar_pipeline).pack(pady=15, ipady=5)
        
        f_log = ttk.LabelFrame(self, text="Monitoring Console")
        f_log.pack(fill='both', expand=True, padx=10, pady=5)
        self.txt_log = tk.Text(f_log, height=14, bg='black', fg='#00FF00', font=('Consolas', 10))
        self.txt_log.pack(fill='both', expand=True, padx=5, pady=5)

    def log(self, msg):
        self.txt_log.insert(tk.END, msg + "\n")
        self.txt_log.see(tk.END)
        self.update_idletasks()

    def AbrirArquivo(self, var):
        p = filedialog.askopenfilename()
        if p: var.set(p)

    def AbrirDiretorio(self, var):
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
        self.log(" -> Loading TRANSIT file (.genes.txt)...")
        df_teste = pd.read_csv(self.var_transit_teste.get(), sep="\t", skiprows=4)
        
        col_gene = 'gene' if 'gene' in df_teste.columns else ('#ORF' if '#ORF' in df_teste.columns else df_teste.columns[0])
        df_teste.rename(columns={col_gene: 'ID_GENE'}, inplace=True)
    
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
            df_m.to_csv(os.path.join(out_dir, "Essentiality_Filter_Control.txt"), sep="\t", index=False)
            return dict(zip(df_m['ID_GENE'], df_m['call_condicao'].fillna("--")))
        else:
            return dict(zip(df_teste['ID_GENE'], df_teste['call'].fillna("--")))

    def extrair_rede_para_csv(self, arquivo_eggnog, cache_rxn, cache_comp, caminho_csv, organism_label):
        self.log(f"-> Processing structured cache for {organism_label}...")
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
            self.log(f"   + Searching {len(novas_rxn)} structural data in KEGG...")
            with ThreadPoolExecutor(max_workers=4) as ex:
                res = list(ex.map(BucaReacoesKEGG, novas_rxn))
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
                            "Substrates": ";".join(subs) if subs else "None",
                            "Enzime": enzima,
                            "Products": ";".join(prods) if prods else "None"
                        })
                        
        df_rede = pd.DataFrame(linhas_rede)
        df_rede.to_csv(caminho_csv, index=False)
        
        novos_comp = [c for c in compostos_detectados if c not in cache_comp or "HierarquiaBRITE" not in cache_comp[c] or "RotaMetabolica" not in cache_comp[c]]
        if novos_comp:
            self.log(f"Updating structural metadata of {len(novos_comp)} compounds via KEGG API...")
            res_lotes = BuscarCompostosKeggLote(novos_comp)
            for c_id, d in res_lotes.items(): cache_comp[c_id] = d
            for c_id in novos_comp:
                if c_id not in cache_comp: 
                    cache_comp[c_id] = {"names": c_id, "RotaMetabolica": [], "HierarquiaBRITE": []}
                
        return df_rede, cache_rxn, cache_comp

    # ==========================================
    # 3. Construction of the association and essentiality network
    # ==========================================

    def EexecutarClassificacaoDaRede(self, df_bac, mapa_essencialidade, cache_comp, out_dir):
        self.log("\n[BUILDING A DECISION ARCHITECTURE FOR METABOLITES]")
        
        todos_mets_brutos = []
        for _, r in df_bac.iterrows():
            if r['Substrates'] != "None": todos_mets_brutos.extend(r['Substrates'].split(';'))
            if r['Products'] != "None": todos_mets_brutos.extend(r['Products'].split(';'))
        
        contagem_hubs = Counter(todos_mets_brutos)
        limite_corte = self.var_corte_hubs.get()
        hubs_removidos = {m for m, c in contagem_hubs.items() if c > limite_corte}
        self.log(f"-> Hubs Filter: {len(hubs_removidos)} overrepresented metabolites excluded (> {limite_corte} connections).")

        enzimas_da_rede = set(df_bac['Enzime'].unique())
        linhas_genes_csv = []
        for enz in enzimas_da_rede:
            classe_ess = mapa_essencialidade.get(enz, "--")
            linhas_genes_csv.append({"Gene": enz, "Essentiality": classe_ess})
        df_genes_mapeados = pd.DataFrame(linhas_genes_csv)
        df_genes_mapeados.to_csv(os.path.join(out_dir, "Gene_Essentiality_Mapping.csv"), index=False)
        mapa_genes_final = dict(zip(df_genes_mapeados['Gene'], df_genes_mapeados['Essentiality']))

        set_substratos = set()
        set_produtos = set()
        for _, row in df_bac.iterrows():
            if row['Substrates'] != "None":
                set_substratos.update([s for s in row['Substrates'].split(';') if s not in hubs_removidos])
            if row['Products'] != "None":
                set_produtos.update([p for p in row['Products'].split(';') if p not in hubs_removidos])
                
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
                enzimas_produtoras = df_bac[df_bac['Products'].str.contains(met, na=False)]['Enzime'].unique()
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
                enzimas_consumidoras = df_bac[df_bac['Substrates'].str.contains(met, na=False)]['Enzime'].unique()
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
                "Topological_Classification": class_topologica,
                "Essentiality_Producers": class_produtores,
                "Essentiality_Consumers": class_consumidores
            })

        df_master_metabolitos = pd.DataFrame(linhas_metabolitos_master)
        df_master_metabolitos.to_csv(os.path.join(out_dir, "Classificacao_Metabolitos_Master.csv"), index=False)
        return df_master_metabolitos, hubs_removidos, mapa_genes_final

    # ==========================================
    # 4. Exporting Results
    # ==========================================

    def exportar_enriquecimento_grupo(self, lista_compounds, cache_comp, prefixo, out_dir):
        """Generation of reports on names, BRITE hierarchies, and PATHWAY."""
        pasta_grupo = os.path.join(out_dir, "Subgroup_Stratification")
        os.makedirs(pasta_grupo, exist_ok=True)
        
        linhas_comp = []
        for c_id in lista_compounds:
            info = cache_comp.get(c_id, {})
            linhas_comp.append({"Metabolito": c_id, "Nome": info.get("names", c_id)})
        df_lista = pd.DataFrame(linhas_comp)
        df_lista.to_csv(os.path.join(pasta_grupo, f"Compound_List_{prefixo}.csv"), index=False, encoding='utf-8')
        
        if not lista_compounds: return
        
        level_counts = {}
        pathway_counts = Counter()
        
        for c_id in lista_compounds:
            info = cache_comp.get(c_id, {})
            for pw in info.get("RotaMetabolica", []):
                pathway_counts[pw] += 1
                
            for path in info.get("HierarquiaBRITE", []):
                for level_idx, categoria in enumerate(path):
                    if level_idx not in level_counts:
                        level_counts[level_idx] = Counter()
                    level_counts[level_idx][categoria] += 1
                    
        linhas_brite = []
        for lvl in sorted(level_counts.keys()):
            for cat, count in level_counts[lvl].most_common():
                linhas_brite.append({
                    "Hierarchy_Level": f"Level{lvl}",
                    "BRITE_Catgory": cat,
                    "Occurrences": count
                })
        df_brite = pd.DataFrame(linhas_brite)
        df_brite.to_csv(os.path.join(pasta_grupo, f"BRITE_Enrichment_{prefixo}.csv"), index=False, encoding='utf-8')
        
        linhas_pw = []
        for pw, count in pathway_counts.most_common():
            linhas_pw.append({"Metabolic_Pathway": pw, "Occurrences": count})
        df_pw = pd.DataFrame(linhas_pw)
        df_pw.to_csv(os.path.join(pasta_grupo, f"PATHWAY_Enrichment_{prefixo}.csv"), index=False, encoding='utf-8')

    def processar_gargalos_finais(self, df_master_met, df_bac, hubs_removidos, mapa_genes, produtos_hospedeiro, cache_comp, out_dir):
        self.log("-> Isolating and classifying biological bottleneck ...")
        
        linhas_gargalos = []
        df_filtrado_gargalo = df_master_met[
            (df_master_met['Essentiality_Consumers'].isin(['ES', 'GD'])) &
            (df_master_met['Essentiality_Producers'].isin(['NE', 'GA', 'SS', '--']))
        ]

        for _, row in df_filtrado_gargalo.iterrows():
            met = row['Metabolito']
            class_prod = row['Essentiality_Producers']
            
            enz_consumidoras = df_bac[df_bac['Substrates'].str.contains(met, na=False)]['Enzime'].unique()
            enz_vulneraveis = [e for e in enz_consumidoras if mapa_genes.get(e, "--") in ['ES', 'GD']]
            
            if class_prod == "SS":
                tipo_gargalo = "External Base (Non-synthesized substrate)"
            else:
                enz_produtoras = df_bac[df_bac['Products'].str.contains(met, na=False)]['Enzime'].unique()
                precursores_vulneraveis = False
                for ep in enz_produtoras:
                    subs_da_ep = df_bac[df_bac['Enzime'] == ep]['Substrates'].values
                    for ListaB in subs_da_ep:
                        if ListaB == "None": continue
                        for s in ListaB.split(';'):
                            if s in hubs_removidos: continue
                            if df_master_met[df_master_met['Metabolito'] == s]['Essentiality_Producers'].isin(['ES', 'GD']).any():
                                precursores_vulneraveis = True
                                
                tipo_gargalo = "Squeezed (Important -> Weak -> Important)" if precursores_vulneraveis else "Internal Base (Weak -> Weak -> Important)"

            esta_na_planta = met in produtos_hospedeiro
            linhas_gargalos.append({
                "Metabolite_Interface": met,
                "Metabolite_Name": cache_comp.get(met, {}).get('names', met),
                "Topological_Classification": row['Topological_Classification'],
                "Bottleneck_Type": tipo_gargalo,
                "Present_In_The_Host": "YES" if esta_na_planta else "NO",
                "Important_Affected_Enzymes": ";".join(enz_vulneraveis)
            })

        df_gargalos_total = pd.DataFrame(linhas_gargalos)
        df_gargalos_total.to_csv(os.path.join(out_dir, "Metabolic_Bottleneck_Mapping.csv"), index=False, encoding='utf-8')

        gargalos_validos = df_gargalos_total[df_gargalos_total['Topological_Classification'].isin(['S', 'M'])]
        
        subgrupos_matriz = {
            "General_ Bottlenecks_SM_Present_in_Host": gargalos_validos[gargalos_validos['Present_In_The_Host'] == "YES"],
            "General_ Bottlenecks_SM_Absent_in_Host": gargalos_validos[gargalos_validos['Present_In_The_Host'] == "NO"],
            
            "S_Present_in_Host": gargalos_validos[(gargalos_validos['Topological_Classification'] == "S") & (gargalos_validos['Present_In_The_Host'] == "YES")],
            "S_Absent_in_Host": gargalos_validos[(gargalos_validos['Topological_Classification'] == "S") & (gargalos_validos['Present_In_The_Host'] == "NO")],
            
            "M_Present_in_Host": gargalos_validos[(gargalos_validos['Topological_Classification'] == "M") & (gargalos_validos['Present_In_The_Host'] == "YES")],
            "M_Absent_in_Host": gargalos_validos[(gargalos_validos['Topological_Classification'] == "M") & (gargalos_validos['Present_In_The_Host'] == "NO")]
        }

        self.log("\n[PERFORMING SUBGROUP ENRICHMENT]")
        relatorio_linhas = [
            "=========================================================",
            " Metabolic Stratification Report  ",
            "=========================================================",
            f"Total of Bottlenecks Analyzed (S+M): {len(gargalos_validos)}\n"
        ]

        for nome_grupo, df_sub in subgrupos_matriz.items():
            IDsLista = df_sub['Metabolite_Interface'].tolist()
            self.exportar_enriquecimento_grupo(IDsLista, cache_comp, nome_grupo, out_dir)
            msg = f"  -> {nome_grupo}: {len(IDsLista)} mapped compounds."
            self.log(msg)
            relatorio_linhas.append(msg)

        with open(os.path.join(out_dir, "Statistical_Summary_Bottlenecks.txt"), "w", encoding='utf-8') as f:
            f.write("\n".join(relatorio_linhas))

    # ==========================================
    # 5. Pipeline Central Execution
    # ==========================================

    def rodar_pipeline(self):
        if not self.var_projeto_dir.get() or not self.var_transit_teste.get() or not self.var_eggnog_bac.get():
            messagebox.showerror("Error", "Fill in the required inputs before running.")
            return
            
        out_dir = os.path.join(self.var_projeto_dir.get(), "Pipeline_Outputs")
        os.makedirs(out_dir, exist_ok=True)
        
        self.txt_log.delete('1.0', tk.END)
        self.log("[STARTING PIPELINE]")

        rxn_cache, comp_cache = self.carregar_banco_master()
        mapa_essencialidade = self.carregar_dados_essencialidade(out_dir)

        df_bac_rede, rxn_cache, comp_cache = self.extrair_rede_para_csv(
            self.var_eggnog_bac.get(), rxn_cache, comp_cache, 
            os.path.join(out_dir, "Bacterial_Cache_Network.csv"), "Bacteria"
        )
        
        produtos_hospedeiro = set()
        if self.var_eggnog_hospedeiro.get():
            df_hosp_rede, rxn_cache, comp_cache = self.extrair_rede_para_csv(
                self.var_eggnog_hospedeiro.get(), rxn_cache, comp_cache,
                os.path.join(out_dir, "Host_Cache_Network.csv"), "Host"
            )
            for _, r in df_hosp_rede.iterrows():
                if r['Products'] != "None": produtos_hospedeiro.update(r['Products'].split(';'))

        self.salvar_banco_master(rxn_cache, comp_cache)

        df_master_met, hubs_removidos, mapa_genes = self.EexecutarClassificacaoDaRede(
            df_bac_rede, mapa_essencialidade, comp_cache, out_dir
        )

        self.processar_gargalos_finais(
            df_master_met, df_bac_rede, hubs_removidos, 
            mapa_genes, produtos_hospedeiro, comp_cache, out_dir
        )
        
        self.log("\n=== COMPLETED! Check the folder 'Pipeline_Outputs/Subgroup_Stratification' ===")
        messagebox.showinfo("Success", "Analysis Complete! Subgroups successfully exported.")

if __name__ == "__main__":
    app = TnSeqMetabolicApp()
    app.mainloop()