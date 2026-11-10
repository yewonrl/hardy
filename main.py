
import streamlit as st
import pandas as pd
import numpy as np
import math
import requests

# ==========================================
# 1. 페이지 설정
# ==========================================
st.set_page_config(
    page_title="열성 유전병 보인자 빈도 예측기",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 하디-바인베르크 법칙 기반 열성 유전병 보인자 빈도 분석기")

st.markdown("""
> 열성 유전병은 질병 관련 열성 대립유전자를 두 개 모두 보유할 때 발병한다.
> 본 프로그램은 하디-바인베르크 법칙을 활용하여 집단 내 정상인,
> 보인자, 환자의 빈도를 추정하고 근친교배가 유전자형 분포에
> 미치는 영향을 분석한다.
""")

st.divider()


# ==========================================
# 2. NCBI API 연동
# ==========================================
try:
    NCBI_API_KEY = st.secrets.get("NCBI_API_KEY", None)
except Exception:
    NCBI_API_KEY = None


def fetch_gene_data(gene_id, api_key=None):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    params = {
        "db": "gene",
        "id": str(gene_id),
        "retmode": "json"
    }

    if api_key:
        params["api_key"] = api_key

    try:
        response = requests.get(
            url,
            params=params,
            timeout=10
        )

        response.raise_for_status()

        data = response.json()
        gene_info = data.get("result", {}).get(
            str(gene_id), {}
        )

        if not gene_info:
            return None

        return {
            "symbol": gene_info.get("name", "Unknown"),
            "description": gene_info.get(
                "description", "정보 없음"
            ),
            "summary": gene_info.get(
                "summary", "설명이 없습니다."
            )
        }

    except requests.exceptions.RequestException:
        return None

    except (ValueError, KeyError):
        return None


# ==========================================
# 3. 유전병 기본 데이터
# ==========================================
preset_data = {
    "낭성 섬유증 (Cystic Fibrosis)": {
        "symbol": "CFTR",
        "gene_id": "1080",
        "default_q2": 0.0004
    },

    "겸상 적혈구 빈혈증 (Sickle Cell)": {
        "symbol": "HBB",
        "gene_id": "3043",
        "default_q2": 0.0025
    },

    "페닐케톤뇨증 (PKU)": {
        "symbol": "PAH",
        "gene_id": "5053",
        "default_q2": 0.0001
    }
}


# ==========================================
# 4. 사이드바 입력 설정
# ==========================================
st.sidebar.header("🔑 1. 유전자 및 API 설정")

disease = st.sidebar.selectbox(
    "분석할 열성 유전병 선택",
    ["직접 입력"] + list(preset_data.keys())
)

default_prevalence = 0.0004
gene_id = None
gene_symbol = "사용자 설정"

if disease != "직접 입력":

    disease_info = preset_data[disease]

    gene_id = disease_info["gene_id"]
    gene_symbol = disease_info["symbol"]

    default_prevalence = disease_info["default_q2"]

    st.sidebar.write(
        f"**연관 유전자:** {gene_symbol}"
    )

    st.sidebar.write(
        f"**NCBI Gene ID:** {gene_id}"
    )

    use_api = st.sidebar.checkbox(
        "NCBI 유전자 정보 조회",
        value=True
    )

    if use_api:

        with st.sidebar.spinner("NCBI 데이터 조회 중..."):

            gene_result = fetch_gene_data(
                gene_id,
                NCBI_API_KEY
            )

        if gene_result:

            st.sidebar.success("NCBI 연동 성공")

            st.sidebar.write(
                f"**유전자 기호:** {gene_result['symbol']}"
            )

            with st.sidebar.expander("유전자 설명"):
                st.write(gene_result["description"])
                st.write(gene_result["summary"])

        else:
            st.sidebar.warning(
                "API 조회에 실패했습니다. 기본 데이터로 진행합니다."
            )

else:

    gene_symbol = st.sidebar.text_input(
        "유전자 기호 입력",
        value="사용자 설정"
    )

    default_prevalence = 0.0004


st.sidebar.divider()

# ==========================================
# 5. 발병률 입력
# ==========================================
st.sidebar.header("⚙️ 2. 발병률 입력")

input_mode = st.sidebar.radio(
    "데이터 입력 방식",
    [
        "기본 프리셋 사용",
        "백분율 (%)",
        "비율 (N명 중 1명)",
        "실제 환자 수 직접 입력"
    ]
)

prevalence = default_prevalence

if input_mode == "백분율 (%)":

    percent = st.sidebar.number_input(
        "환자 발생 백분율 (%)",
        min_value=0.000001,
        max_value=100.0,
        value=max(default_prevalence * 100, 0.000001),
        step=0.01
    )

    prevalence = percent / 100

elif input_mode == "비율 (N명 중 1명)":

    ratio = st.sidebar.number_input(
        "몇 명 중 1명꼴인가요?",
        min_value=1,
        value=max(
            int(1 / default_prevalence), 1
        ),
        step=100
    )

    prevalence = 1 / ratio

elif input_mode == "실제 환자 수 직접 입력":

    sample_pop = st.sidebar.number_input(
        "조사 대상 인구수",
        min_value=1,
        value=100000,
        step=10000
    )

    sample_patients = st.sidebar.number_input(
        "실제 환자 수",
        min_value=0,
        max_value=int(sample_pop),
        value=min(
            int(sample_pop * default_prevalence),
            int(sample_pop)
        ),
        step=1
    )

    prevalence = sample_patients / sample_pop


st.sidebar.divider()


# ==========================================
# 6. 인구 규모 및 근친교배 계수
# ==========================================
st.sidebar.header("🔬 3. 인구 규모 및 심화 변수")

total_population = st.sidebar.number_input(
    "분석 대상 집단 총인구수",
    min_value=1000,
    max_value=100000000,
    value=100000,
    step=10000
)

use_inbreeding = st.sidebar.checkbox(
    "근친교배 계수(F) 적용"
)

if use_inbreeding:

    f_coef = st.sidebar.slider(
        "근친교배 계수(F)",
        min_value=0.0,
        max_value=1.0,
        value=0.0625,
        step=0.01
    )

else:
    f_coef = 0.0


# ==========================================
# 7. 하디-바인베르크 계산 함수
# ==========================================
def calculate_genotype(prevalence, F):

    if not 0 <= prevalence <= 1:
        raise ValueError("발병률은 0~1 사이여야 합니다.")

    if not 0 <= F <= 1:
        raise ValueError("F는 0~1 사이여야 합니다.")

    # 발병률 = q² + Fpq
    # p = 1-q를 대입하여 q 계산

    if F == 1:
        q = prevalence

    else:
        q = (
            -F + math.sqrt(
                F**2 + 4 * (1-F) * prevalence
            )
        ) / (2 * (1-F))

    p = 1 - q

    freq_AA = p**2 + F*p*q
    freq_Aa = 2*p*q*(1-F)
    freq_aa = q**2 + F*p*q

    return {
        "p": p,
        "q": q,
        "AA": freq_AA,
        "Aa": freq_Aa,
        "aa": freq_aa
    }


# ==========================================
# 8. 인구수 계산 함수
# ==========================================
def calculate_population(population, result):

    frequencies = [
        result["AA"],
        result["Aa"],
        result["aa"]
    ]

    raw_counts = [
        population * x for x in frequencies
    ]

    counts = [
        math.floor(x) for x in raw_counts
    ]

    remaining = population - sum(counts)

    order = sorted(
        range(3),
        key=lambda i: raw_counts[i] - counts[i],
        reverse=True
    )

    for i in order[:remaining]:
        counts[i] += 1

    return counts


# ==========================================
# 9. 계산 실행
# ==========================================
if not 0 <= prevalence <= 1:

    st.error("발병률은 0~100% 사이로 입력해야 합니다.")
    st.stop()

result = calculate_genotype(
    prevalence,
    f_coef
)

p = result["p"]
q = result["q"]

freq_AA = result["AA"]
freq_Aa = result["Aa"]
freq_aa = result["aa"]

num_normal, num_carriers, num_affected = (
    calculate_population(total_population, result)
)

ratio_carrier_to_patient = (
    freq_Aa / freq_aa if freq_aa > 0 else 0
)


# ==========================================
# 10. 기초 대립유전자 빈도
# ==========================================
st.subheader("📊 1. 기초 대립유전자 빈도")

col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "열성 대립유전자 빈도(q)",
    f"{q:.5f}"
)

col2.metric(
    "정상 대립유전자 빈도(p)",
    f"{p:.5f}"
)

col3.metric(
    "근친교배 계수(F)",
    f"{f_coef:.4f}"
)

col4.metric(
    "환자 1명당 보인자 비율",
    f"{ratio_carrier_to_patient:.1f}명"
)

st.divider()


# ==========================================
# 11. 예상 유전자형 구성
# ==========================================
st.subheader(
    f"👥 2. 인구 {total_population:,}명 기준 예상 구성"
)

col1, col2, col3 = st.columns(3)

col1.metric(
    "정상인 (AA)",
    f"{num_normal:,}명",
    f"{freq_AA*100:.4f}%"
)

col2.metric(
    "보인자 (Aa)",
    f"{num_carriers:,}명",
    f"{freq_Aa*100:.4f}%"
)

col3.metric(
    "열성 환자 (aa)",
    f"{num_affected:,}명",
    f"{freq_aa*100:.4f}%"
)

st.divider()


# ==========================================
# 12. 대립유전자 빈도 변화 그래프
# ==========================================
st.subheader(
    "📈 3. 열성 대립유전자 빈도에 따른 유전자형 변화"
)

q_range = np.linspace(0.001, 0.5, 100)
p_range = 1 - q_range

df_graph = pd.DataFrame({

    "열성 대립유전자 빈도(q)": q_range,

    "정상인(AA)": (
        p_range**2 + p_range*q_range*f_coef
    ),

    "보인자(Aa)": (
        2*p_range*q_range*(1-f_coef)
    ),

    "환자(aa)": (
        q_range**2 + p_range*q_range*f_coef
    )
})

st.line_chart(
    df_graph.set_index("열성 대립유전자 빈도(q)")
)

st.divider()


# ==========================================
# 13. 근친교배 계수 비교
# ==========================================
st.subheader("🧬 4. 근친교배 계수별 유전자형 비교")

F_values = [0.0, 0.0625, 0.125, 0.25, 0.5, 1.0]

comparison_data = []

for F in F_values:

    r = calculate_genotype(prevalence, F)

    comparison_data.append({

        "근친교배 계수(F)": F,

        "정상인 비율(%)": r["AA"] * 100,

        "보인자 비율(%)": r["Aa"] * 100,

        "환자 비율(%)": r["aa"] * 100

    })

df_comparison = pd.DataFrame(comparison_data)

st.dataframe(
    df_comparison,
    use_container_width=True
)

st.bar_chart(
    df_comparison.set_index("근친교배 계수(F)")[
        ["정상인 비율(%)", "보인자 비율(%)", "환자 비율(%)"]
    ]
)

st.divider()


# ==========================================
# 14. 생명과학적 결과 해석
# ==========================================
st.subheader("🧪 5. 생명과학 및 의학적 분석")

st.info(
    f"""
    **학술적 해석**

    설정한 조건에서 열성 유전병 환자(aa)는
    전체 집단의 **{freq_aa*100:.4f}% ({num_affected:,}명)**로 추정됩니다.

    반면 특별한 증상이 나타나지 않는 보인자(Aa)는
    전체 집단의 **{freq_Aa*100:.4f}% ({num_carriers:,}명)**로 추정됩니다.

    따라서 열성 유전병은 실제 환자 수만으로
    집단 내 열성 대립유전자의 분포를 파악하기 어려우며,
    보인자 빈도를 함께 분석하는 것이 중요합니다.
    """
)

if f_coef > 0:

    st.warning(
        f"""
        **근친교배 영향 분석**

        근친교배 계수 F={f_coef:.4f}를 적용한 결과,
        이형접합자인 보인자(Aa)의 비율은
        {(freq_Aa*100):.4f}%이며,
        열성 동형접합자인 환자(aa)의 비율은
        {(freq_aa*100):.4f}%입니다.

        근친교배는 이형접합자 비율을 감소시키고
        동형접합자 비율을 증가시키므로
        열성 유전병의 발현 가능성에 영향을 줍니다.
        """
    )

else:

    st.success(
        """
        **무작위 교배 집단 분석**

        근친교배 계수를 0으로 설정하여
        하디-바인베르크 평형의 기본 조건에 따른
        유전자형 분포를 계산하였습니다.
        """
    )

st.divider()


# ==========================================
# 15. 결과 데이터 표 및 CSV 다운로드
# ==========================================
st.subheader("📋 6. 세부 분석 결과")

df_summary = pd.DataFrame({

    "구분": [
        "정상인 (AA)",
        "보인자 (Aa)",
        "열성 환자 (aa)"
    ],

    "유전자형 수식": [
        "p² + Fpq",
        "2pq(1-F)",
        "q² + Fpq"
    ],

    "집단 내 비율 (%)": [
        f"{freq_AA*100:.4f}%",
        f"{freq_Aa*100:.4f}%",
        f"{freq_aa*100:.4f}%"
    ],

    "예상 인구수 (명)": [
        f"{num_normal:,}",
        f"{num_carriers:,}",
        f"{num_affected:,}"
    ]

})

st.dataframe(
    df_summary,
    use_container_width=True
)

csv_data = df_summary.to_csv(
    index=False
).encode("utf-8-sig")

st.download_button(
    label="📥 분석 결과 CSV 다운로드",
    data=csv_data,
    file_name="genetic_carrier_analysis.csv",
    mime="text/csv"
)

st.divider()

st.caption(
    "본 프로그램은 하디-바인베르크 법칙에 기반한 "
    "집단 수준의 이론적 빈도 추정 시뮬레이션입니다."
)
