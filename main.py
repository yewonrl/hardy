import streamlit as st
import pandas as pd
import numpy as np

# Page Setup
st.set_page_config(
    page_title="하디-바인베르크 & FST 시뮬레이터",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 무작위 유전적 부동과 고정지수(F_ST) 분화 시뮬레이터")
st.markdown("하디-바인베르크 평형 조건이 깨지는 **유전적 부동(Genetic Drift)과 집단 격리** 환경에서 $F_{ST}$ 지수 변화를 탐구합니다.")
st.divider()

# Sidebar Setup
st.sidebar.header("⚙️ 집단유전학 변수 설정")

p_init = st.sidebar.slider("초기 대립유전자 A 빈도 (p)", min_value=0.1, max_value=0.9, value=0.5, step=0.05)
num_subpops = st.sidebar.slider("격리된 하위 집단 개수 (Subpopulations, k)", min_value=2, max_value=10, value=5, step=1)
pop_size = st.sidebar.slider("하위 집단별 개체 수 (N)", min_value=10, max_value=1000, value=50, step=10)
generations = st.sidebar.number_input("추적할 세대 수 (Generations)", min_value=10, max_value=300, value=100, step=10)

st.sidebar.divider()
st.sidebar.subheader("🔀 이주율 (Gene Flow)")
migration_rate = st.sidebar.slider("집단 간 이주율 (m)", min_value=0.00, max_value=0.20, value=0.01, step=0.01)

# Simulation Logic
np.random.seed(42)  # 재현성을 위한 시드 고정

p_matrix = np.zeros((generations + 1, num_subpops))
p_matrix[0, :] = p_init

fst_list = []
hs_list = []
ht_list = []

# Generation Loop
for g in range(generations):
    current_p = p_matrix[g, :]
    
    # 1. HWE 연산 및 FST 산출
    mean_p = np.mean(current_p)
    H_S = np.mean(2 * current_p * (1 - current_p))
    H_T = 2 * mean_p * (1 - mean_p)
    
    F_ST = (H_T - H_S) / H_T if H_T > 0 else 0.0
    
    hs_list.append(round(H_S, 4))
    ht_list.append(round(H_T, 4))
    fst_list.append(round(max(0.0, F_ST), 4))
    
    # 2. 이주(Gene Flow) 효과 반영
    p_migrated = (1 - migration_rate) * current_p + migration_rate * mean_p
    
    # 3. 무작위 유전적 부동 (이항 샘플링)
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

# Output UI
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
