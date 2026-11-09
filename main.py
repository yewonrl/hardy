import math

def print_intro():
    print("=" * 70)
    print(" 🧬 하디-바인베르크 법칙: 열성 유전병 보인자 빈도 예측 프로그램 🧬")
    print("-" * 70)
    print("[참고 지문]")
    print("유전 질환의 경우 특정 유전자가 두 개 모두 있을 때 질병이 나타나는")
    print("열성 유전병에 대해, 해당 유전자를 하나만 지닌 보인자가 인구 중")
    print("얼마나 되는지를 추정할 수 있습니다. 이는 질병의 실제 환자 빈도와 함께 계산되어,")
    print("해당 유전자가 집단에 얼마나 퍼져 있는지를 예측하는 데 도움이 됩니다.")
    print("=" * 70)

def get_detailed_input():
    while True:
        print("\n[상세 조건 입력 설정]")
        print("어떤 방식으로 환자 데이터를 입력하시겠습니까?")
        print("1. 백분율(%)로 입력 (예: 0.01% -> 0.01 입력)")
        print("2. 비율로 입력 (예: 10,000명 중 1명 -> 10000 입력)")
        print("3. 전체 인구수와 실제 환자 수로 직접 입력")
        
        choice = input("\n입력 방식 선택 (1, 2, 3): ").strip()
        
        try:
            q_squared = 0
            population = 0
            
            if choice == '1':
                percent = float(input("👉 환자 발생 백분율(%)을 입력하세요: "))
                q_squared = percent / 100
                population = int(input("👉 결과를 추정할 전체 집단 인구수를 입력하세요(예: 100000): "))
                break
                
            elif choice == '2':
                ratio = int(input("👉 몇 명 중 1명 꼴로 환자가 발생하나요? (숫자만 입력): "))
                q_squared = 1 / ratio
                population = int(input("👉 결과를 추정할 전체 집단 인구수를 입력하세요(예: 100000): "))
                break
                
            elif choice == '3':
                population = int(input("👉 전체 집단의 인구수를 입력하세요: "))
                patients = int(input(f"👉 {population}명 중 실제 환자 수는 몇 명입니까?: "))
                q_squared = patients / population
                break
                
            else:
                print("⚠️ 잘못된 입력입니다. 1, 2, 3 중에서 선택해주세요.")
                
        except ValueError:
            print("⚠️ 숫자만 정확히 입력해주세요. 다시 시도합니다.")
            
    return q_squared, population

def calculate_hwe(q_squared, population):
    # 하디-바인베르크 방정식 계산
    # p = 정상 대립유전자 빈도, q = 열성(질환) 대립유전자 빈도
    
    q = math.sqrt(q_squared)
    p = 1 - q
    
    # 유전자형 빈도 계산
    p_squared = p ** 2       # 정상(우성 동형접합) 빈도
    two_pq = 2 * p * q       # 보인자(이형접합) 빈도
    
    # 실제 인구수에 따른 사람 수 추산 (소수점 반올림)
    normal_count = round(p_squared * population)
    carrier_count = round(two_pq * population)
    patient_count = round(q_squared * population)
    
    # 결과 출력
    print("\n" + "=" * 70)
    print(" 📊 분석 결과 리포트")
    print("=" * 70)
    
    print(f"▶ [대립유전자 빈도]")
    print(f" - 정상 유전자(A) 빈도 (p): {p:.4f} ({p*100:.2f}%)")
    print(f" - 질환 유전자(a) 빈도 (q): {q:.4f} ({q*100:.2f}%)")
    print(f"   (💡 질환 유전자가 집단 내에 약 {q*100:.2f}% 퍼져 있음을 의미합니다.)")
    
    print(f"\n▶ [유전자형 빈도 및 예상 인구 분포 (총 {population:,}명 기준)]")
    print(f" - 환자 (aa, {q_squared*100:.4f}%): 약 {patient_count:,}명")
    print(f" - 보인자 (Aa, {two_pq*100:.2f}%): 약 {carrier_count:,}명")
    print(f" - 정상 (AA, {p_squared*100:.2f}%): 약 {normal_count:,}명")
    print("-" * 70)
    
    # 지문의 핵심 내용을 반영한 해석
    carrier_to_patient_ratio = carrier_count / patient_count if patient_count > 0 else float('inf')
    print("💡 [분석 결과 해석]")
    print(f"실제 겉으로 드러나는 환자는 {patient_count:,}명에 불과하지만,")
    print(f"질병 유전자를 숨기고 있는 보인자는 {carrier_count:,}명으로 환자보다 약 {carrier_to_patient_ratio:.0f}배 많습니다.")
    print("이는 희귀 유전 질환일수록 보인자의 형태로 집단 내에 유전자가 광범위하게")
    print("퍼져 숨어있음을 하디-바인베르크 법칙을 통해 수학적으로 입증하는 결과입니다.")
    print("=" * 70)

# 프로그램 실행
if __name__ == "__main__":
    print_intro()
    q_squared, population = get_detailed_input()
    calculate_hwe(q_squared, population)
