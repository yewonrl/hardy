import streamlit as st
import pandas as pd
import numpy as np
import math
import requests

# Page Setup
st.set_page_config(
    page_title="하디-바인베르크 유전병 보인자 추정기",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 하디-바인베르크 법칙 기반 열성 유전병 보인자 빈도 분석기")
st.markdown(
    "> *\"유전 질환의 경우 특정 유전자가 두 개 모두 있을 때 질병이 나타나는 열성 유전병에 대해, "
    "해당 유전자를 하나만 지닌 보인자가 인구 중 얼마나 되는 지를 추정할 수 있다. "
    "이는 질병의 실제 환자 빈도와 함께 계산되어, 해당 유전자가 집단에 얼마나 퍼져 있는 지를 예측하는 데 도움이 된다.\"*"
)
st.divider()

# ---------------------------------------------------------
# Secrets & API Handling (NCBI Entrez API)
# ---------------------------------------------------------
try:
    ncbi_api_key = st.secrets["NCBI_API_KEY"]
    has_secret_key = True
except Exception:
    ncbi_api_key = None
    has_secret_key = False

def fetch_gene_info(gene_name, api_key=None):
    """NCBI Entrez API를 통해 열성 유전병 관련 유전자 ID 조회"""
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "gene",
            "term": f"{gene_name}[Gene Name] AND human[Organism]",
            "retmode": "json"
        }
        if api_key:
            params["api_key"] = api_key
            
        res = requests.get(url, params=params, timeout=4)
        if res.status_code == 200:
            data = res.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            if id_list:
                return True, id_list[0], f"NCBI Gene ID: {id_list[0]} 연동 성공"
        return False, None, "NCBI 데이터베이스에서 해당 유전자를 찾을 수 없습니다."
    except Exception:
        return False, None, "네트워크 통신 예외 발생 (기본 모드로 동작합니다)."

# ---------------------------------------------------------
# Sidebar Setup
# ---------------------------------------------------------
st.sidebar.header("🔑 1. API 보안 및 유전자 조회")

if has_secret_key:
    st.sidebar.success("🔒 Streamlit Secrets: NCBI_API_KEY 연동 완료")
else:
    st.sidebar.info("💡 Secrets 미설정: Public API 모드로 동작합니다.")

preset_disease = st.sidebar.selectbox(
    "주요 열성 유전병 예시 선택",
    [
        "직접 입력 (사용자 설정)",
        "낭성 섬유증 (Cystic Fibrosis, CFTR 유전자)",
        "페닐케톤뇨증 (PKU, PAH 유전자)",
        "겸상 적혈구 빈혈증 (Sickle Cell Anemia, HBB 유전자)",
        "테이-삭스병 (Tay-Sachs Disease, HEXA 유전자)"
    ]
)

# Preset data (환자 발병 빈도 q^2)
preset_data = {
    "낭성 섬유증 (Cystic Fibrosis, CFTR 유전자)": ("CFTR", 0.0004),      # 1/2500
    "페닐케톤뇨증 (PKU, PAH 유전자)": ("PAH", 0.0001),                    # 1/10000
    "겸상 적혈구 빈혈증 (Sickle Cell Anemia, HBB 유전자)": ("HBB", 0.0025), # 1/400
    "테이-삭스병 (Tay-Sachs Disease, HEXA 유전자)": ("HEXA", 0.00028)    # 1/3600
}

if preset_disease != "직접 입력 (사용자 설정)":
    gene_symbol, default_q2 = preset_data[preset_disease]
    q2_input = default_q2
    st.sidebar.write(f"**연관 유전자:** `{gene_symbol}`")
else:
    gene_symbol = st.sidebar.text_input("유전자 기호 입력", value="CFTR")
    q2_input = st.sidebar.number_input("열성 환자 발병률 (q²)", min_value=0.000001, max_value=0.200000, value=0.000400, step=0.000100, format="%.6f")

use_api = st.sidebar.checkbox("NCBI API 유전자 정보 파싱 연동", value=True)
if use_api and gene_symbol:
    api_success, g_id, api_msg = fetch_gene_info(gene_symbol, ncbi_api_key)
    if api_success:
        st.sidebar.success(api_msg)
    else:
        st.sidebar.warning(api_msg)

st.sidebar.divider()
st.sidebar.header("⚙️ 2. 집단 규모 및 조건 설정")
total_population = st.sidebar.number_input("분석 대상 집단 인구수 (명)", min_value=1000, max_value=100000000, value=100000, step=10000)

# ---------------------------------------------------------
# Hardy-Weinberg Calculation Logic
# ---------------------------------------------------------
# q^2 = 환자 빈도 (aa)
# q = sqrt(q^2) -> 열성 대립유전자 빈도
# p = 1 - q    -> 정상 대립유전자 빈도
# p^2 = 정상인 (AA)
# 2pq = 보인자 (Aa)

q = math.sqrt(q2_input)
p = 1.0 - q

p2_freq = p ** 2
two_pq_freq = 2 * p * q

num_affected = int(round(total_population * q2_input))
num_carriers = int(round(total_population * two_pq_freq))
num_normal = total_population - num_affected - num_carriers

# Carrier to Patient Ratio
ratio_carrier_to_patient = round(two_pq_freq / q2_input, 1) if q2_input > 0 else 0

# ---------------------------------------------------------
# UI Main Output
# ---------------------------------------------------------
st.subheader("📊 1. 유전자 빈도 및 인구 집단 추정 메트릭")

col1, col2, col3, col4 = st.columns(4)
col1.metric("열성 대립유전자 빈도 (q)", f"{q:.4f}")
col2.metric("정상 대립유전자 빈도 (p)", f"{p:.4f}")
col3.metric("추정 보인자 비율 (2pq)", f"{two_pq_freq * 100:.2f}%")
col4.metric("환자 1명당 보인자 비율", f"{ratio_carrier_to_patient:.1f} 배")

st.divider()

st.subheader(f"👥 2. 인구 {total_population:,}명 기준 예상 유전자형 구성")

col_a, col_b, col_c = st.columns(3)
col_a.metric("정상인 (AA)", f"{num_normal:,} 명", f"{p2_freq*100:.2f}%")
col_b.metric("보인자 (Aa, 증상 없음)", f"{num_carriers:,} 명", f"{two_pq_freq*100:.2f}%")
col_c.metric("열성 환자 (aa)", f"{num_affected:,} 명", f"{q2_input*100:.4f}%")

st.divider()

st.subheader("📈 3. 열성 대립유전자 빈도(q) 변화에 따른 유전자형 비율 곡선")

# Graph Simulation Range
q_range = np.linspace(0.001, 0.5, 100)
p_range = 1.0 - q_range

df_graph = pd.DataFrame({
    "열성 대립유전자 빈도 (q)": q_range,
    "정상인 비율 (p²)": p_range ** 2,
    "보인자 비율 (2pq)": 2 * p_range * q_range,
    "환자 비율 (q²)": q_range ** 2
})

chart_df = df_graph.set_index("열성 대립유전자 빈도 (q)")
st.line_chart(chart_df)

st.divider()

st.subheader("🧪 4. 세특용 생명과학 및 의학적 해석 보고서")

st.info(
    f"💡 **학술 해석:**  \n"
    f"선택한 조건에서 실제 **열성 질환 환자(aa)**는 집단의 **{q2_input*100:.4f}%({num_affected:,}명)**에 불과하지만, "
    f"증상이 전혀 없는 **잠재적 보인자(Aa)**는 전체 인구의 **{two_pq_freq*100:.2f}%({num_carriers:,}명)**에 달합니다.  \n"
    f"즉, **환자 1명이 발생할 때 집단 내에는 약 {ratio_carrier_to_patient:.1f}명의 보인자가 숨어 있음**을 의미합니다."
)

st.warning(
    "🚨 **보인자 빈도 추정의 중요성:**  \n"
    "희귀 열성 유전병이라 하더라도 대립유전자 빈도 수식($q = \\sqrt{\\text{발병률}}$)의 특성상, "
    "열성 유전자가 집단 내에서 완전히 사라지지 않고 상당수 보인자(Aa)의 형태로 보존됩니다. "
    "이는 유전 상담, 산전 유전자 검사 및 유전병 예방 프로그램 수립에 있어 결정적인 수학적 근거가 됩니다."
)

st.divider()
st.subheader("📋 5. 세부 수치 데이터 시뮬레이션 표")

df_summary = pd.DataFrame({
    "구분": ["정상인 (AA)", "보인자 (Aa)", "환자 (aa)"],
    "유전자형 수식": ["p²", "2pq", "q²"],
    "집단 내 비율 (%)": [f"{p2_freq*100:.4f}%", f"{two_pq_freq*100:.4f}%", f"{q2_input*100:.6f}%"],
    "추정 인구수 (명)": [f"{num_normal:,}", f"{num_carriers:,}", f"{num_affected:,}"]
})

st.dataframe(df_summary, use_container_width=True)
