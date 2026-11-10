
import streamlit as st
import pandas as pd
import numpy as np
import requests
import math


# =========================================================
# 1. 페이지 설정
# =========================================================

st.set_page_config(
    page_title="하디-바인베르크 유전병 보인자 예측기",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 하디-바인베르크 법칙 기반 열성 유전병 보인자 빈도 분석기")

st.markdown("""
> 열성 유전병은 특정 열성 대립유전자를 두 개 모두 물려받았을 때
> 발병하는 유전 질환이다. 반면 하나의 열성 대립유전자만 가진
> 보인자는 일반적으로 증상이 나타나지 않는다.
>
> 본 프로그램은 하디-바인베르크 법칙을 기반으로 집단 내
> 열성 유전병 환자와 잠재적 보인자의 빈도를 추정하고,
> 근친교배가 유전자형 분포에 미치는 영향을 분석한다.
""")

st.divider()


# =========================================================
# 2. NCBI API 설정
# =========================================================

try:
    NCBI_API_KEY = st.secrets["NCBI_API_KEY"]
except Exception:
    NCBI_API_KEY = None


@st.cache_data(ttl=3600)
def fetch_gene_data(gene_id, api_key=None):
    """NCBI Gene DB에서 유전자 정보 조회"""

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

        gene_info = data.get(
            "result", {}
        ).get(str(gene_id), {})

        if not gene_info or "error" in gene_info:
            return None

        return {
            "symbol": gene_info.get("name", "Unknown"),
            "description": gene_info.get("description", ""),
            "summary": gene_info.get("summary", "")
        }

    except (requests.RequestException, ValueError, KeyError):
        return None


# =========================================================
# 3. 유전병 프리셋 데이터
# =========================================================

PRESET_DATA = {
    "낭성 섬유증 (Cystic Fibrosis)": {
        "symbol": "CFTR",
        "gene_id": "1080",
        "rs_id": "rs113993960",
        "default_q2": 0.0004
    },

    "겸상 적혈구 빈혈증 (Sickle Cell)": {
        "symbol": "HBB",
        "gene_id": "3043",
        "rs_id": "rs334",
        "default_q2": 0.0025
    },

    "페닐케톤뇨증 (PKU)": {
        "symbol": "PAH",
        "gene_id": "5053",
        "rs_id": "rs5030858",
        "default_q2": 0.0001
    }
}


# =========================================================
# 4. 핵심 계산 함수
# =========================================================

def calculate_genotype(prevalence, F):
    """
    발병률과 근친교배 계수 F를 이용하여
    대립유전자 및 유전자형 빈도를 계산한다.

    발병률 = q² + Fpq
    p + q = 1
    """

    if not 0 <= prevalence <= 1:
        raise ValueError("발병률은 0~1 사이여야 합니다.")

    if not 0 <= F <= 1:
        raise ValueError("근친교배 계수는 0~1 사이여야 합니다.")

    if F == 1:
        q = prevalence

    else:
        q = (
            -F + math.sqrt(
                F**2 + 4 * (1 - F) * prevalence
            )
        ) / (2 * (1 - F))

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


def allocate_population(result, population):
    """반올림 오차를 보정하여 총인구수를 유지한다."""

    frequencies = np.array([
        result["AA"],
        result["Aa"],
        result["aa"]
    ])

    raw_counts = frequencies * population

    counts = np.floor(raw_counts).astype(int)

    remainder = population - counts.sum()

    order = np.argsort(
        -(raw_counts - counts)
    )

    for i in range(remainder):
        counts[order[i]] += 1

    return counts


def make_comparison_data(q, population, f_values):
    """동일한 q에서 F값별 유전자형 구성 비교"""

    p = 1 - q
    data = []

    for F in f_values:

        result = {
            "AA": p**2 + F*p*q,
            "Aa": 2*p*q*(1-F),
            "aa": q**2 + F*p*q
        }

        counts = allocate_population(
            result,
            population
        )

        data.append({
            "근친교배 계수(F)": F,
            "정상인 비율(%)": result["AA"] * 100,
            "보인자 비율(%)": result["Aa"] * 100,
            "환자 비율(%)": result["aa"] * 100,
            "정상인 수(명)": counts[0],
            "보인자 수(명)": counts[1],
            "환자 수(명)": counts[2]
        })

    return pd.DataFrame(data)


# =========================================================
# 5. 사이드바 - 유전자 선택 및 API
# =========================================================

st.sidebar.header("🔑 1. 유전자 및 API 설정")

if NCBI_API_KEY:
    st.sidebar.success("NCBI API Key 연동 완료")
else:
    st.sidebar.info("NCBI Public API 모드")

preset_disease = st.sidebar.selectbox(
    "분석할 열성 유전병 선택",
    ["직접 입력 (사용자 설정)"] + list(PRESET_DATA.keys())
)

default_q2 = 0.0004
gene_symbol = ""
gene_id = ""
rs_id = ""

if preset_disease != "직접 입력 (사용자 설정)":

    disease_info = PRESET_DATA[preset_disease]

    gene_symbol = disease_info["symbol"]
    gene_id = disease_info["gene_id"]
    rs_id = disease_info["rs_id"]
    default_q2 = disease_info["default_q2"]

    st.sidebar.markdown(
        f"**연관 유전자:** `{gene_symbol}`"
    )

    st.sidebar.markdown(
        f"**NCBI Gene ID:** `{gene_id}`"
    )

    st.sidebar.markdown(
        f"**관련 변이 ID:** `{rs_id}`"
    )

    use_api = st.sidebar.checkbox(
        "NCBI 유전자 정보 조회",
        value=True
    )

    if use_api:

        with st.sidebar.spinner("NCBI 데이터 조회 중..."):

            api_result = fetch_gene_data(
                gene_id,
                NCBI_API_KEY
            )

        if api_result:

            st.sidebar.success(
                f"유전자 조회 완료: {api_result['symbol']}"
            )

            with st.sidebar.expander("유전자 상세 정보"):

                st.write(api_result["description"])
                st.write(api_result["summary"])

        else:
            st.sidebar.warning(
                "NCBI 조회 실패. 기본 데이터로 진행합니다."
            )

else:

    gene_symbol = st.sidebar.text_input(
        "유전자 기호 입력",
        value="사용자 입력 유전자"
    )

st.sidebar.caption(
    "NCBI의 유전자 정보는 조회하지만, SNP의 MAF를 "
    "질병 대립유전자 빈도로 간주하지는 않습니다."
)

st.sidebar.divider()


# =========================================================
# 6. 사이드바 - 발병 빈도 입력
# =========================================================

st.sidebar.header("⚙️ 2. 환자 발생 빈도 설정")

st.sidebar.caption(
    "프리셋 수치는 실제 임상 통계가 아닌 시뮬레이션용 기본값입니다."
)

input_mode = st.sidebar.radio(
    "데이터 입력 방식",
    [
        "프리셋 기본값 사용",
        "백분율 (%)",
        "비율 (N명 중 1명)",
        "실제 환자 수 직접 입력"
    ]
)

q2_input = default_q2

if input_mode == "백분율 (%)":

    percent = st.sidebar.number_input(
        "환자 발생 백분율 (%)",
        min_value=0.0,
        max_value=100.0,
        value=float(default_q2 * 100),
        step=0.01,
        format="%.6f"
    )

    q2_input = percent / 100

elif input_mode == "비율 (N명 중 1명)":

    default_ratio = max(1, round(1 / default_q2))

    ratio = st.sidebar.number_input(
        "몇 명 중 1명꼴인가요?",
        min_value=1,
        value=default_ratio,
        step=100
    )

    q2_input = 1 / ratio

elif input_mode == "실제 환자 수 직접 입력":

    sample_pop = st.sidebar.number_input(
        "조사 대상 인구수 (명)",
        min_value=1,
        value=100000,
        step=10000
    )

    sample_patients = st.sidebar.number_input(
        "실제 환자 수 (명)",
        min_value=0,
        max_value=int(sample_pop),
        value=min(
            int(sample_pop * default_q2),
            int(sample_pop)
        ),
        step=1
    )

    q2_input = sample_patients / sample_pop

st.sidebar.divider()


# =========================================================
# 7. 사이드바 - 인구 규모 및 근친교배 변수
# =========================================================

st.sidebar.header("🔬 3. 심화 분석 변수")

total_population = st.sidebar.number_input(
    "분석 대상 집단의 총인구수 (명)",
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
        step=0.01,
        help="F=0은 무작위 교배를 가정한 상태입니다."
    )

else:
    f_coef = 0.0

st.sidebar.divider()


# =========================================================
# 8. 전체 계산 실행
# =========================================================

if not 0 <= q2_input <= 1:

    st.error("발병률은 0~100% 범위여야 합니다.")
    st.stop()

result = calculate_genotype(
    q2_input,
    f_coef
)

p = result["p"]
q = result["q"]

freq_AA = result["AA"]
freq_Aa = result["Aa"]
freq_aa = result["aa"]

num_normal, num_carriers, num_affected = allocate_population(
    result,
    int(total_population)
)

ratio_carrier_to_patient = (
    freq_Aa / freq_aa
    if freq_aa > 0
    else None
)


# =========================================================
# 9. 기초 대립유전자 빈도 분석
# =========================================================

st.subheader("📊 1. 기초 대립유전자 빈도 분석")

col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "열성 대립유전자 빈도(q)",
    f"{q:.6f}"
)

col2.metric(
    "정상 대립유전자 빈도(p)",
    f"{p:.6f}"
)

col3.metric(
    "근친교배 계수(F)",
    f"{f_coef:.4f}"
)

col4.metric(
    "환자 1명당 보인자 비율",
    f"{ratio_carrier_to_patient:.1f}배"
    if ratio_carrier_to_patient is not None
    else "계산 불가"
)

st.divider()


# =========================================================
# 10. 유전자형별 예상 인구수
# =========================================================

st.subheader(
    f"👥 2. 전체 인구 {total_population:,}명 기준 유전자형 구성"
)

col_a, col_b, col_c = st.columns(3)

col_a.metric(
    "정상인 (AA)",
    f"{num_normal:,}명",
    f"{freq_AA*100:.4f}%"
)

col_b.metric(
    "보인자 (Aa)",
    f"{num_carriers:,}명",
    f"{freq_Aa*100:.4f}%"
)

col_c.metric(
    "열성 유전병 환자 (aa)",
    f"{num_affected:,}명",
    f"{freq_aa*100:.4f}%"
)

st.info(
    f"전체 집단에서 예상되는 보인자는 {num_carriers:,}명이며, "
    f"전체 인구의 {freq_Aa*100:.2f}%를 차지합니다."
)

st.divider()


# =========================================================
# 11. 대립유전자 빈도에 따른 유전자형 그래프
# =========================================================

st.subheader(
    "📈 3. 열성 대립유전자 빈도(q)에 따른 유전자형 변화"
)

q_range = np.linspace(0, 0.5, 200)
p_range = 1 - q_range

df_graph = pd.DataFrame({

    "정상인 (AA)": (
        p_range**2 + f_coef*p_range*q_range
    ) * 100,

    "보인자 (Aa)": (
        2*p_range*q_range*(1-f_coef)
    ) * 100,

    "환자 (aa)": (
        q_range**2 + f_coef*p_range*q_range
    ) * 100

}, index=q_range)

df_graph.index.name = "열성 대립유전자 빈도(q)"

st.line_chart(df_graph)

st.caption(
    f"현재 근친교배 계수 F = {f_coef:.4f} 기준"
)

st.divider()


# =========================================================
# 12. 근친교배 계수별 비교
# =========================================================

st.subheader(
    "🧬 4. 근친교배 계수(F)에 따른 유전자형 변화"
)

st.markdown("""
동일한 열성 대립유전자 빈도(q)를 유지한 상태에서
근친교배 계수만 변화시켜 유전자형 분포를 비교한다.
""")

f_values = [0, 0.0625, 0.125, 0.25, 0.5, 1.0]

df_F = make_comparison_data(
    q,
    int(total_population),
    f_values
)

st.dataframe(
    df_F.style.format({
        "근친교배 계수(F)": "{:.4f}",
        "정상인 비율(%)": "{:.4f}",
        "보인자 비율(%)": "{:.4f}",
        "환자 비율(%)": "{:.4f}"
    }),
    use_container_width=True,
    hide_index=True
)

st.markdown("#### 보인자 및 환자 비율 비교")

chart_F = df_F.set_index("근친교배 계수(F)")[
    ["보인자 비율(%)", "환자 비율(%)"]
]

st.line_chart(chart_F)

st.info("""
**분석 해석**

동일한 q에서 근친교배 계수 F가 증가하면
이형접합자인 보인자(Aa)의 비율은 감소하고,
열성 동형접합자인 환자(aa)의 비율은 증가한다.
""")

st.divider()


# =========================================================
# 13. 근친교배 적용 전후 그래프 비교
# =========================================================

st.subheader(
    "📈 5. 근친교배 적용 전후 유전자형 변화 비교"
)

q_compare = np.linspace(0.001, 0.5, 200)
p_compare = 1 - q_compare

df_compare = pd.DataFrame({

    "F=0 보인자": (
        2*p_compare*q_compare
    ) * 100,

    "현재 F 보인자": (
        2*p_compare*q_compare*(1-f_coef)
    ) * 100,

    "F=0 환자": (
        q_compare**2
    ) * 100,

    "현재 F 환자": (
        q_compare**2 + f_coef*p_compare*q_compare
    ) * 100

}, index=q_compare)

df_compare.index.name = "열성 대립유전자 빈도(q)"

st.markdown("#### 보인자 비율 비교")

st.line_chart(
    df_compare[
        ["F=0 보인자", "현재 F 보인자"]
    ]
)

st.markdown("#### 열성 유전병 환자 비율 비교")

st.line_chart(
    df_compare[
        ["F=0 환자", "현재 F 환자"]
    ]
)

st.caption(
    f"현재 적용된 근친교배 계수: F = {f_coef:.4f}"
)

st.divider()


# =========================================================
# 14. 생명과학 및 의학적 해석
# =========================================================

st.subheader(
    "🧪 6. 생명과학 및 의학적 해석 보고서"
)

st.info(
    f"""
**[학술적 분석 결과]**

설정한 발병률과 근친교배 계수에 따라 분석한 결과,
열성 유전병 환자(aa)는 집단의
**{freq_aa*100:.4f}%({num_affected:,}명)**로 추정되었다.

반면 열성 대립유전자를 하나만 가진 보인자(Aa)는
전체의 **{freq_Aa*100:.2f}%({num_carriers:,}명)**로 나타났다.

이는 실제 환자 수만으로는 파악하기 어려운 잠재적 보인자
집단의 규모를 수학적으로 추정할 수 있음을 보여준다.
"""
)

if f_coef > 0:

    st.warning(
        f"""
**[근친교배 효과 분석]**

현재 근친교배 계수 F = {f_coef:.4f}가 적용되었다.

근친교배가 증가하면 이형접합자(Aa)의 빈도는 감소하고,
동형접합자(AA, aa)의 빈도는 증가한다.

따라서 동일한 대립유전자 빈도에서도
열성 유전병이 발현되는 동형접합자(aa)의 비율이
높아질 수 있다.
"""
    )

else:

    st.success(
        """
**[보인자 빈도 추정의 의의]**

열성 유전병은 환자 수가 적더라도 집단 내에
해당 열성 대립유전자가 보인자 형태로 존재할 수 있다.

따라서 환자 빈도뿐 아니라 보인자 빈도를 함께 분석하는 것은
집단의 유전적 구성과 유전 질환의 분포를 이해하는 데 중요하다.
"""
    )

st.divider()


# =========================================================
# 15. 유전자형별 세부 데이터 표
# =========================================================

st.subheader("📋 7. 유전자형별 세부 분석 데이터")

df_summary = pd.DataFrame({

    "유전자형": [
        "정상인 (AA)",
        "보인자 (Aa)",
        "환자 (aa)"
    ],

    "계산식": [
        "p² + Fpq",
        "2pq(1-F)",
        "q² + Fpq"
    ],

    "유전자형 빈도": [
        freq_AA,
        freq_Aa,
        freq_aa
    ],

    "집단 내 비율(%)": [
        freq_AA*100,
        freq_Aa*100,
        freq_aa*100
    ],

    "예상 인구수(명)": [
        num_normal,
        num_carriers,
        num_affected
    ]

})

st.dataframe(
    df_summary.style.format({
        "유전자형 빈도": "{:.6f}",
        "집단 내 비율(%)": "{:.4f}"
    }),
    use_container_width=True,
    hide_index=True
)

st.divider()


# =========================================================
# 16. CSV 다운로드
# =========================================================

st.subheader("📥 8. 분석 결과 다운로드")

df_export = df_summary.copy()

df_export["근친교배 계수(F)"] = f_coef
df_export["열성 대립유전자 빈도(q)"] = q
df_export["정상 대립유전자 빈도(p)"] = p
df_export["분석 대상 인구수"] = int(total_population)

csv_data = df_export.to_csv(
    index=False,
    encoding="utf-8-sig"
)

st.download_button(
    label="📥 분석 결과 CSV 다운로드",
    data=csv_data,
    file_name="hardy_weinberg_analysis.csv",
    mime="text/csv",
    use_container_width=True
)

st.markdown("#### 근친교배 비교 데이터 다운로드")

comparison_csv = df_F.to_csv(
    index=False,
    encoding="utf-8-sig"
)

st.download_button(
    label="📥 근친교배 비교 CSV 다운로드",
    data=comparison_csv,
    file_name="inbreeding_comparison.csv",
    mime="text/csv",
    use_container_width=True
)

st.divider()


# =========================================================
# 17. 수학적 원리 및 분석의 한계
# =========================================================

st.subheader("📚 9. 분석 원리 및 한계")

with st.expander("하디-바인베르크 법칙과 근친교배 계산 원리"):

    st.markdown(r"""
    ### 기본 하디-바인베르크 법칙

    \[
    p+q=1
    \]

    \[
    p^2+2pq+q^2=1
    \]

    ### 근친교배 계수 적용

    \[
    AA=p^2+Fpq
    \]

    \[
    Aa=2pq(1-F)
    \]

    \[
    aa=q^2+Fpq
    \]

    ### 발병률을 이용한 q 역산

    \[
    r=q^2+Fq(1-q)
    \]

    \[
    q=\frac{-F+\sqrt{F^2+4(1-F)r}}{2(1-F)}
    \]

    단, F=1이면 q=r로 계산한다.
    """)

with st.expander("프로그램의 분석 한계 및 활용 분야"):

    st.markdown("""
    ### 분석의 한계

    1. 본 프로그램은 단일 유전자 좌위에 대한 이론적 모델이다.
    2. 실제 집단에서는 자연선택, 돌연변이, 유전자 이동 및 집단 구조의 영향을 받는다.
    3. 실제 유전병의 발병률은 지역과 인구집단에 따라 달라질 수 있다.
    4. 프리셋 발병률은 교육용 시뮬레이션 기본값이며 실제 임상 통계가 아니다.
    5. SNP의 MAF는 질병 대립유전자 빈도와 반드시 일치하지 않는다.
    6. 근친교배 계수는 집단의 유전적 구조를 단순화한 값이며,
       모든 개인의 실제 근친 정도를 의미하지 않는다.
    7. 실제 유전병은 여러 변이와 유전적 배경의 영향을 받을 수 있다.

    ### 활용 분야

    - 집단유전학
    - 유전 질환 연구
    - 생명정보학
    - 유전자 검사 및 유전 상담의 기초 이해
    - 생명공학 연구
    """)

st.caption(
    "본 프로그램은 교육 및 탐구 목적의 이론적 시뮬레이션 도구이며, "
    "개인의 질병 위험이나 보인자 여부를 진단하지 않습니다."
)
