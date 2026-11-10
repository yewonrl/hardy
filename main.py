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
# Secrets & API Handling (NCBI 다중 DB 연동)
# ---------------------------------------------------------
try:
    ncbi_api_key = st.secrets["NCBI_API_KEY"]
    has_secret_key = True
except Exception:
    ncbi_api_key = None
    has_secret_key = False

# [새로 추가된 통합 API 함수]
def fetch_ncbi_combined_data(gene_id, rs_id, api_key=None):
    """NCBI의 gene DB와 snp DB를 모두 호출하여 정보와 대립유전자 빈도를 가져오는 함수"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    result_data = {"success": False, "symbol": "", "description": "", "q_value": None}
    
    # 1. db=gene 호출 (유전자 정보)
    try:
        gene_params = {"db": "gene", "id": gene_id, "retmode": "json"}
        if api_key: gene_params["api_key"] = api_key
        gene_resp = requests.get(base_url, params=gene_params, timeout=4)
        if gene_resp.status_code == 200:
            gene_info = gene_resp.json().get("result", {}).get(str(gene_id), {})
            result_data["symbol"] = gene_info.get("name", "Unknown")
            result_data["description"] = gene_info.get("summary", "설명이 없습니다.")
            result_data["success"] = True
    except Exception:
        pass # 실패 시 기본값 유지

    # 2. db=snp 호출 (대립유전자 빈도)
    if rs_id:
        snp_id_only = rs_id.replace("rs", "")
        try:
            snp_params = {"db": "snp", "id": snp_id_only, "retmode": "json"}
            if api_key: snp_params["api_key"] = api_key
            snp_resp = requests.get(base_url, params=snp_params, timeout=4)
            if snp_resp.status_code == 200:
                snp_info = snp_resp.json().get("result", {}).get(snp_id_only, {})
                if "global_maf" in snp_info:
                    maf_str = snp_info["global_maf"]
                    if "=" in maf_str and "/" in maf_str:
                        result_data["q_value"] = float(maf_str.split("=")[1].split("/")[0])
        except Exception:
            pass

    return result_data

# ---------------------------------------------------------
# Sidebar Setup: 1. API & 유전자 선택
# ---------------------------------------------------------
st.sidebar.header("🔑 1. API 보안 및 유전자 조회")

if has_secret_key:
    st.sidebar.success("🔒 Streamlit Secrets 연동 완료")
else:
    st.sidebar.info("💡 Secrets 미설정: Public API 모드")

# [수정된 프리셋 데이터: Gene ID와 rsID 추가]
preset_data = {
    "낭성 섬유증 (Cystic Fibrosis)": {"symbol": "CFTR", "gene_id": "1080", "rs_id": "rs113993960", "default_q2": 0.0004},
    "겸상 적혈구 빈혈증 (Sickle Cell)": {"symbol": "HBB", "gene_id": "3043", "rs_id": "rs334", "default_q2": 0.0025},
    "페닐케톤뇨증 (PKU)": {"symbol": "PAH", "gene_id": "5053", "rs_id": "rs5030858", "default_q2": 0.0001}
}

preset_disease = st.sidebar.selectbox(
    "주요 열성 유전병 예시 선택",
    ["직접 입력 (사용자 설정)"] + list(preset_data.keys())
)

default_q2 = 0.0004
api_q_value = None

if preset_disease != "직접 입력 (사용자 설정)":
    disease_info = preset_data[preset_disease]
    st.sidebar.write(f"**연관 유전자:** `{disease_info['symbol']}` (ID: {disease_info['gene_id']})")
    
    use_api = st.sidebar.checkbox("NCBI API 유전자 정보 및 발병률 파싱", value=True)
    
    if use_api:
        with st.sidebar.spinner('Gene & SNP 연동 중...'):
            api_result = fetch_ncbi_combined_data(disease_info['gene_id'], disease_info['rs_id'], ncbi_api_key)
            
        if api_result["success"]:
            st.sidebar.success(f"NCBI 정보 연동 완료 ({api_result['symbol']})")
            api_q_value = api_result["q_value"]
            
            if api_q_value is not None and api_q_value > 0:
                st.sidebar.info(f"🌐 실시간 글로벌 빈도(q): {api_q_value:.4f}")
                default_q2 = api_q_value ** 2 # API로 가져온 q값을 q^2으로 역산하여 기본값 대체
            else:
                st.sidebar.warning("SNP 빈도 데이터 없음. 프리셋 수치를 사용합니다.")
                default_q2 = disease_info['default_q2']
        else:
            st.sidebar.error("API 통신 실패. 기본 모드로 동작합니다.")
            default_q2 = disease_info['default_q2']
    else:
        default_q2 = disease_info['default_q2']
else:
    st.sidebar.text_input("유전자 기호 입력", value="직접 입력 모드")
    default_q2 = 0.0004

st.sidebar.divider()

# ---------------------------------------------------------
# Sidebar Setup: 2. 디테일한 환자 빈도 입력 방식
# ---------------------------------------------------------
st.sidebar.header("⚙️ 2. 환자 발생 빈도 상세 입력")
input_mode = st.sidebar.radio(
    "데이터 입력 방식을 선택하세요:",
    ["프리셋/API 연동값 사용", "백분율 (%)", "비율 (N명 중 1명)", "실제 환자 수 직접 입력"]
)

q2_input = default_q2

if input_mode == "백분율 (%)":
    percent = st.sidebar.number_input("환자 발생 백분율(%)", min_value=0.000001, max_value=100.0, value=default_q2*100, step=0.01, format="%.4f")
    q2_input = percent / 100.0
elif input_mode == "비율 (N명 중 1명)":
    ratio = st.sidebar.number_input("몇 명 중 1명 꼴로 발생하나요? (명)", min_value=1, value=int(1/default_q2) if default_q2 > 0 else 10000, step=100)
    q2_input = 1.0 / ratio
elif input_mode == "실제 환자 수 직접 입력":
    sample_pop = st.sidebar.number_input("조사 대상 인구수 (명)", min_value=1, value=100000, step=10000)
    sample_patients = st.sidebar.number_input("실제 환자 수 (명)", min_value=0, max_value=sample_pop, value=int(100000*default_q2), step=1)
    q2_input = sample_patients / sample_pop if sample_pop > 0 else 0

st.sidebar.divider()

# ---------------------------------------------------------
# Sidebar Setup: 3. 인구 규모 및 심화 환경 변수
# ---------------------------------------------------------
st.sidebar.header("🔬 3. 인구 규모 및 심화 변수")
total_population = st.sidebar.number_input("결과를 추산할 집단 총 인구수 (명)", min_value=1000, max_value=100000000, value=100000, step=10000)

use_inbreeding = st.sidebar.checkbox("심화: 근친교배 계수(F) 적용")
if use_inbreeding:
    f_coef = st.sidebar.slider(
        "근친교배 계수 (F)", 0.0, 1.0, 0.0625, 0.01, 
        help="0: 완전 무작위 교배, 0.0625: 사촌 간 교배, 0.25: 남매 간 교배"
    )
else:
    f_coef = 0.0

# ---------------------------------------------------------
# Hardy-Weinberg Calculation Logic (F 계수 포함)
# ---------------------------------------------------------
q = math.sqrt(q2_input)
p = 1.0 - q

freq_AA = (p ** 2) + (p * q * f_coef)
freq_Aa = (2 * p * q) * (1 - f_coef)
freq_aa = (q ** 2) + (p * q * f_coef)

num_affected = int(round(total_population * freq_aa))
num_carriers = int(round(total_population * freq_Aa))
num_normal = total_population - num_affected - num_carriers

ratio_carrier_to_patient = round(freq_Aa / freq_aa, 1) if freq_aa > 0 else 0

# ---------------------------------------------------------
# UI Main Output
# ---------------------------------------------------------
st.subheader("📊 1. 기초 대립유전자 빈도 메트릭")

col1, col2, col3, col4 = st.columns(4)
col1.metric("열성 대립유전자 빈도 (q)", f"{q:.4f}")
col2.metric("정상 대립유전자 빈도 (p)", f"{p:.4f}")
col3.metric("근친교배 계수 (F)", f"{f_coef:.4f}" if use_inbreeding else "미적용")
col4.metric("환자 1명당 보인자 비율", f"{ratio_carrier_to_patient:.1f} 배")

st.divider()

st.subheader(f"👥 2. 인구 {total_population:,}명 기준 예상 유전자형 구성")

col_a, col_b, col_c = st.columns(3)
col_a.metric("정상인 (AA)", f"{num_normal:,} 명", f"{freq_AA*100:.2f}%")
col_b.metric("보인자 (Aa, 증상 없음)", f"{num_carriers:,} 명", f"{freq_Aa*100:.2f}%")
col_c.metric("열성 환자 (aa)", f"{num_affected:,} 명", f"{freq_aa*100:.4f}%")

st.divider()

st.subheader("📈 3. 열성 대립유전자 빈도(q) 변화에 따른 유전자형 비율 곡선")

q_range = np.linspace(0.001, 0.5, 100)
p_range = 1.0 - q_range

df_graph = pd.DataFrame({
    "열성 대립유전자 빈도 (q)": q_range,
    "정상인 비율 (AA)": (p_range ** 2) + (p_range * q_range * f_coef),
    "보인자 비율 (Aa)": (2 * p_range * q_range) * (1 - f_coef),
    "환자 비율 (aa)": (q_range ** 2) + (p_range * q_range * f_coef)
})

chart_df = df_graph.set_index("열성 대립유전자 빈도 (q)")
st.line_chart(chart_df)

st.divider()

st.subheader("🧪 4. 세특용 생명과학 및 의학적 해석 보고서")

st.info(
    f"💡 **학술 해석:**  \n"
    f"설정된 조건에서 실제 겉으로 드러나는 **열성 질환 환자(aa)**는 집단의 **{freq_aa*100:.4f}%({num_affected:,}명)**에 불과하지만, "
    f"증상이 전혀 없는 **잠재적 보인자(Aa)**는 전체 인구의 **{freq_Aa*100:.2f}%({num_carriers:,}명)**에 달합니다.  \n"
    f"즉, **환자 1명이 발생할 때 집단 내에는 약 {ratio_carrier_to_patient:.1f}명의 보인자가 겉보기엔 정상인처럼 숨어 있음**을 의미합니다."
)

if use_inbreeding and f_coef > 0:
    st.error(
        f"🚨 **근친교배(Inbreeding) 효과 분석:**  \n"
        f"현재 근친교배 계수(F={f_coef})가 적용되어 하디-바인베르크 평형이 수정되었습니다. "
        f"무작위 교배(F=0) 상황에 비해 이형접합자(Aa)의 비율이 줄어들고, 동형접합자(aa)의 비율이 증가하여 "
        f"**열성 유전병이 겉으로 발현될 확률이 수학적으로 높아졌음**을 시뮬레이션을 통해 증명할 수 있습니다."
    )
else:
    st.warning(
        "🚨 **보인자 빈도 추정의 중요성:**  \n"
        "희귀 열성 유전병이라 하더라도 대립유전자 빈도 수식($q = \\sqrt{\\text{발병률}}$)의 특성상, "
        "열성 유전자가 집단 내에서 완전히 사라지지 않고 상당수 보인자(Aa)의 형태로 보존됩니다. "
        "이는 유전 상담, 산전 유전자 검사 및 보건 정책 수립에 있어 결정적인 수학적 근거가 됩니다."
    )

st.divider()
st.subheader("📋 5. 세부 수치 데이터 시뮬레이션 표")

df_summary = pd.DataFrame({
    "구분": ["정상인 (AA)", "보인자 (Aa)", "환자 (aa)"],
    "유전자형 수식 (F적용)": ["p² + pqF", "2pq(1-F)", "q² + pqF"],
    "집단 내 비율 (%)": [f"{freq_AA*100:.4f}%", f"{freq_Aa*100:.4f}%", f"{freq_aa*100:.6f}%"],
    "추정 인구수 (명)": [f"{num_normal:,}", f"{num_carriers:,}", f"{num_affected:,}"]
})

st.dataframe(df_summary, use_container_width=True)
