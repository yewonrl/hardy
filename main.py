import streamlit as st
import pandas as pd
import numpy as np
import requests

# Page Setup
st.set_page_config(
    page_title="하디-바인베르크 & FST 시뮬레이터",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 무작위 유전적 부동과 고정지수(F_ST) 분화 시뮬레이터")
st.markdown("NCBI Open API 데이터 연동 및 인구유전학 수식을 결합하여 **유전적 부동과 집단 격리** 환경을 탐구합니다.")
st.divider()

# ---------------------------------------------------------
# Secrets & API Handling
# ---------------------------------------------------------
# Streamlit Secrets에서 NCBI API 키 불러오기 (없을 시 기본 퍼블릭 키 처리)
try:
    ncbi_api_key = st.secrets["NCBI_API_KEY"]
    has_secret_key = True
except Exception:
    ncbi_api_key = None
    has_secret_key = False

def fetch_ncbi_gene_summary(gene_symbol, api_key=None):
    """NCBI Entrez Direct API를 호출하여 유전자 데이터 조회"""
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "gene",
            "term": f"{gene_symbol}[Gene Name] AND human[Organism]",
            "retmode": "json"
        }
        if api_key:
            params["api_key"] = api_key
            
        response = requests.get(url, params=params, timeout=4)
        if response.status_code == 200:
            data = response.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            if id_list:
                gene_id = id_list[0]
                return True, gene_id, f"NCBI Gene ID: {gene_id} 연동 성공"
            else:
                return False, None, "유전자 검색 결과가 없어 기본 대립유전자 빈도를 사용합니다."
        return False, None, f"API 응답 실패 (Status: {response.status_code})"
    except Exception as e:
        return False, None, f"네트워크 통신 예외 발생 (기본값으로 가상 실행합니다)"

# ---------------------------------------------------------
# Sidebar Setup
# ---------------------------------------------------------
st.sidebar.header("🔑 1. API 보안 및 유전자 검색")

if has_secret_key:
    st.sidebar.success("🔒 Streamlit Secrets: NCBI_API_KEY 인증 성공")
else:
    st.sidebar.info("💡 Secrets 미설정 상태: Public API 모드로 동작합니다.")

use_ncbi_api = st.sidebar.checkbox("NCBI Entrez API 유전자 검색 연동", value=False)

p_init = 0.50

if use_ncbi_api:
    target_gene = st.sidebar.text_input("검색할 인간 유전자 기호 (예: BRCA1, TP53)", value="BRCA1")
    success, g_id, api_msg = fetch_ncbi_gene_summary(target_gene, ncbi_api_key)
    
    if success:
        st.sidebar.success(api_msg)
        # 유전자 ID 해시값을 활용한 임의 대립유전자 초기 빈도 산출 (0.2~0.8)
        p_init = round(0.2 + (int(g_id) % 60) / 100.0, 2)
        st.sidebar.info(f"선택 유전자 기반 대립유전자(A) 초기 빈도(p): **{p_init}**")
    else:
        st.sidebar.warning(api_msg)
        p_init = st.sidebar.slider("초기 대립유전자 A 빈도 (p)", min_value=0.1, max_value=0.9, value=0.5, step=0.05)
else:
    p_init = st.sidebar.slider("초기 대립유전자 A 빈도 (p)", min_value=0.1, max_value=0.9, value=0.5, step=0.05)

st.sidebar.divider()
st.sidebar.header("⚙️ 2. 집단유전학 변수 설정")

num_subpops = st.sidebar.slider("격리된 하위 집단 개수 (k)", min_value=2, max_value=10, value=5, step=1)
pop_size = st.sidebar.slider("하위 집단별 개체 수 (N)", min_value=10, max_value=1000, value=50, step=10)
generations = st.sidebar.number_input("추적할 세대 수 (Generations)", min_value=10, max_value=300, value=100, step=10)
migration_rate = st.sidebar.slider("집단 간 이주율 (m)", min_value=0.00, max_value=0.20, value=0.01, step=0.01)

# ---------------------------------------------------------
# Simulation Logic
# ---------------------------------------------------------
np.random.seed(42)

p_matrix = np.zeros((generations + 1, num_subpops))
p_matrix[0, :] = p_init

fst_list = []
hs_list = []
ht_list = []

for g in range(generations):
    current_p = p_matrix[g, :]
    
    # 1. HWE 및 FST 연산
    mean_p = np.mean(current_p)
    H_S = np.mean(2 * current_p * (1 - current_p))
    H_T = 2 * mean_p * (1 - mean_p)
    
    F_ST = (H_T - H_S) / H_T if H_T > 0 else 0.0
    
    hs_list.append(round(H_S, 4))
    ht_list.append(round(H_T, 4))
    fst_list.append(round(max(0.0, F_ST), 4))
    
    # 2. 이주 및 무작위 유전적 부동 연산
    p_migrated = (1 - migration_rate) * current_p + migration_rate * mean_p
    next_p = np.random.binomial(2 * pop_size, p_migrated) / (2 * pop_size)
    p_matrix[g + 1, :] = next_p

# 마지막 세대 FST 산출
final_p = p_matrix[-1, :]
mean_p_final = np.mean(final_p)
H_S_final = np.mean(2 * final_p * (1 - final_p))
H_T_final = 2 * mean_p_final * (1 - mean_p_final)
F_ST_final = (H_T_final - H_S_final) / H_T_final if H_T_final > 0 else 0.0

hs_list.append(round(H_S_final, 4))
ht_list.append(round(H_T_final, 4))
fst_list.append(round(max(0.0, F_ST_final), 4))

# Dataframe Assembly
df_p = pd.DataFrame(p_matrix, columns=[f"집단 {i+1}" for i in range(num_subpops)])
df_p["세대"] = list(range(generations + 1))

df_stats = pd.DataFrame({
    "세대": list(range(generations + 1)),
    "F_ST (고정지수)": fst_list,
    "H_S (집단내 이형접합도)": hs_list,
    "H_T (전체 이형접합도)": ht_list
})

# ---------------------------------------------------------
# Output UI
# ---------------------------------------------------------
st.subheader("📊 1. 최종 세대 유전적 분화 결과 요약")

col1, col2, col3, col4 = st.columns(4)
col1.metric("초기 p 빈도", f"{p_init:.2f}")
col2.metric("최종 F_ST (고정지수)", f"{fst_list[-1]:.4f}")
col3.metric("최종 H_S (하위집단 이형성)", f"{hs_list[-1]:.4f}")
col4.metric("최종 H_T (전체집단 이형성)", f"{ht_list[-1]:.4f}")

st.divider()

st.subheader("📈 2. 세대별 하위 집단 대립유전자 빈도(p) 및 F_ST 추이")

tab1, tab2 = st.tabs(["하위 집단별 p 빈도 변화", "F_ST (고정지수) 변화 곡선"])

with tab1:
    chart_p = df_p.set_index("세대")
    st.line_chart(chart_p)

with tab2:
    chart_fst = df_stats.set_index("세대")[["F_ST (고정지수)"]]
    st.line_chart(chart_fst)

st.divider()

st.subheader("🧪 3. 인구유전학 학술 메커니즘 분석")

final_fst_val = fst_list[-1]

if final_fst_val < 0.05:
    st.info("💡 **유전적 분화 미미 (F_ST < 0.05):** 활발한 이주(Gene Flow) 또는 큰 집단 크기로 인해 하위 집단 간 유전적 차이가 거의 없으며 하디-바인베르크 평형 상태에 가깝습니다.")
elif final_fst_val < 0.15:
    st.warning("💡 **중등도 유전적 분화 (0.05 ≤ F_ST < 0.15):** 유전적 부동에 의해 하위 집단 간 대립유전자 빈도 차이가 벌어지기 시작했습니다.")
elif final_fst_val < 0.25:
    st.warning("🚨 **높은 유전적 분화 (0.15 ≤ F_ST < 0.25):** 지리적 격리와 부동 효과로 인해 집단 간 유전적 교류가 차단되어 독자적인 유전자 풀을 형성하고 있습니다.")
else:
    st.error("🚨 **극도로 높은 유전적 분화 (F_ST ≥ 0.25):** 대립유전자 고정(Fixation) 현상이 관찰되며, 집단 간 유전적 차이가 심화되어 종 분화(Speciation) 가능성이 존재합니다.")

st.divider()
st.subheader("📋 4. 세대별 상세 수치 데이터")
st.dataframe(df_stats, use_container_width=True)
