# converse-KAM 3D — 위상·동역학·수송 지표의 비교

**Paul, Hudson & Helander (2022)가 "향후 논문"으로 남긴 비교를 완성한다.**

---

## 배경 — 4년 묵은 구멍

**Paul, Hudson & Helander (2022, JPP 88, 905880107) 결론:**

> "our metric appears to be related to the converse KAM approach...
> We expect that the effective volume of parallel diffusion might agree with
> such a calculation in the limit of small perpendicular diffusion.
> **We reserve such a comparison for future publication.**"

**★ 2022년에 적혔고, 4년째 나오지 않았다.**

---

## 무대 — Paul의 네 자기장

```
χ(ψ,θ,ζ) = ι'ψ²/2 + Σ_{m,n} ε_{m,n} · ψ(ψ − ψ̄) · cos(mθ − nζ)
```

**임계 겹침 (Chirikov S = 1.000), separatrix가 ρ = 1/8, 7/8:**

| | 공명 n/m | 공명 수 |
|---|---|---|
| **m=4** | [1/4, 2/4, 3/4] | 3 |
| **m=12** | [2/12, ..., 10/12] | 9 |
| **m=20** | [3/20, ..., 17/20] | 15 |
| **m=36** | [5/36, ..., 31/36] | 27 |

**★ 넷 다 KAM 곡면이 전멸했다. 그런데 열 수송이 완전히 다르다.**

---

## ★★ 최종 결과

**같은 자기장, 같은 격자, 같은 코어(ρ ∈ [0.25, 0.75]), 각 축의 원본 코드 (κ⊥=1e-6):**

| | m=4 | m=12 | m=20 | m=36 | 구분? |
|---|---|---|---|---|---|
| **[위상] converse-KAM 비존재** | 99.7% | 100% | 100% | 100% | **★ 못 함** |
| **[위상] t_c 중앙값 (ζ)** | 18.5 | 20.7 | 19.5 | 20.0 | **★ 못 함** |
| **[동역학] WBA dig 중앙값** | 1.16 | 1.00 | 1.10 | 1.72 | **★ 못 함** |
| **[동역학] 카오스 분율 (dig<5)** | 0.87 | 0.88 | 0.88 | 0.77 | **★ 못 함** |
| **★ [수송] V_PD** | **0.897** | 0.834 | 0.697 | **0.423** | **★ 함** |
| **★ [수송] ΔT** | **0.039** | 0.060 | 0.081 | **0.113** | **★ 함** |

### ★ 읽는 법

**① 위상과 동역학은 넷을 구분하지 못한다.**
- converse-KAM: 코어에서 **100% 포화**. 넷 다 "KAM 없음."
- WBA: dig ≈ 1, 카오스 분율 ~80–88%. 넷 다 카오스.

**② 수송은 구분한다.**
- **V_PD: 0.897 → 0.423** (m=4가 가장 많이 샘)
- **★ ΔT: 0.039 → 0.113** (m=36이 온도 기울기를 가장 잘 유지)

**★ ③ 그리고 결정적으로 — m=36은 공명이 27개다. 가장 많이 깨졌다.**
**★ 그런데 가장 잘 절연한다.**

> **★ "가장 깨진 자기장이 가장 잘 절연한다."**

---

## ★★ 그리고 Paul의 예상 하나가 틀렸다

**Paul (2022) §3의 두 진술:**

| | 진술 | 검정 결과 |
|---|---|---|
| **①** | "converse-KAM은 넷을 구분 못 할 것" | **★ 맞다** (C2) |
| **②** | "κ⊥ → 0에서 converse-KAM과 V_PD가 일치할 것" | **★ 틀리다** (C3b v2) |

**★ 측정 — r(−t_c ↔ V_PD), κ⊥ = 1e-4 / 1e-5 / 1e-6:**

```
  m=4:   0.116 → 0.022 → 0.043    ★ 커지지 않음
  m=36:  0.002 → −0.001 → 0.007   ★ 거의 0
  ★ 네 자기장 전부 소-κ⊥에서 커지지 않는다.
```

### ★ 왜 ②가 틀린가 — 구조적 이유

**★ 두 진술이 양립할 수 없다:**

```
만약 converse-KAM과 V_PD가 일치한다면
→ V_PD가 넷을 구분하므로
→ ★ converse-KAM도 넷을 구분해야 한다
→ ★ 그런데 안 한다 (100% 포화)
```

**★ 그리고 더 근본적으로:**

> **converse-KAM은 코어에서 완전히 포화되어 있다 (비존재 100%).**
> **★ 즉 이진 관측량의 공간 분산이 0이다. 상관을 계산할 대상이 없다.**

**★ "일치할 것"이 없다. 비교할 변량이 없다.**
(소-κ⊥에서의 "일치"는 그냥 양쪽이 포화한 것 — converse-KAM → 100%, V_PD → 1.)

---

## ★ 그래서 — 무엇을 보였나

> **위상 지표(converse-KAM)도 동역학 지표(WBA)도,
> 수송 지표(V_PD)가 보는 것을 놓친다.**
>
> **★ 그리고 그것은 도구의 결함이 아니라,
> 자기장만으로는 수송이 결정되지 않기 때문이다.**

**★ V_PD는 온도장을 푼다. 그리고 온도장은 κ⊥를 통해 플라즈마 상태를 담는다.**
**★ converse-KAM과 WBA는 자기장만 본다. κ⊥를 모른다.**

**★ 이것이 Vlad et al. (2002)의 λ_e가 말하는 것이다:**
> 같은 자기장이라도 플라즈마 상태가 확산계수를 몇 자릿수 바꾼다.

**★ L_c는 분자일 뿐이고, 분모(λ_e)는 플라즈마에서 온다.**

---

## 진행 기록

| | 내용 | 커밋 |
|---|---|---|
| **C0** | KMM·Paul 자기장 구현 + Poincaré 검증 | e268182 |
| **C1** | converse-KAM 3D 구현 + 섬 면적 검증 (0.3%) | c35facb |
| **C2** | converse-KAM을 Paul 자기장 4개에 → **구분 못 함** | 393328e |
| **C3a** | WBA를 같은 자기장에 → **구분 못 함** | bf5ca06 |
| **C3b** | V_PD 1차 (★ 자기장 오류로 SUPERSEDED) | 08df601 |
| **A** | 정합성 확인 → **★ 자기장 모델 불일치 발견** | 823c942 |
| **★ C3b v2** | Paul 명세 자기장으로 재실행 + 세 축 비교 완성 | 8a52ca1 |

상세: `records/RECORD_C0..C3a`, `records/RECORD_C3b_v2_vpd.md`(현행),
`records/RECORD_C3b_vpd.md`(SUPERSEDED), `records/RECORD_A*_*.md`, `records/future_questions.md`.

---

## ★★ 배운 함정

### ① 격자 이산화 (비단조 오차) — C2

섬 면적의 격자 오차가 **비단조 ±0.5%** 요동.
**★ 진단: converse-KAM을 아예 쓰지 않고 해석적 참 섬만 셀 카운트해도
똑같은 요동** → 100% 이산화.

> **★ 격자 오차가 단조 감소하지 않으면, 수렴이 아니라 요동이다.**

### ② 유한시간 착각 — C2, C3a에서 두 번

**C2:** 3000 회전에서 m=36의 카오스가 파편화된 것처럼 보임
→ 10⁵ 회전에서 전 범위를 돔. **cantori 병목.**

**★ C3a:** T=1000에서 m=36의 dig가 가장 높음
→ **★ 그리고 그것이 정확히 우리가 원하던 답(cantori 가설)이었다**
→ T=5000에서 씻겨나감. **곡선이 교차.**

> **★ 유한시간에서 본 것이 무한시간의 진실이 아니다.**
> **★ 그리고 원하던 답이 유한시간에서 나오면 더 의심하라.**

### ③ 눈금 부재 (게이트) — C3a

**★ dig ≈ 1이 "낮다"고 말하려면 "높다"가 뭔지 알아야 한다.**
게이트(단일 공명에서 섬 내부 dig ≈ 12–14)가 그 눈금이다.

> **★ 대비군 없이 측정값을 해석하지 말라.**

### ★★ ④ 검증 맹점 — C3b (가장 교활했다)

**C3b 1차가 Paul 명세를 안 따랐다** (ψ(ψ−ψ̄) 포락선 누락, 상수 진폭 사용).

**★ 그런데 검증을 전부 통과했다** (ε_crit=√κ⊥/2 한 자릿수 이내, V_PD ∝ √ε, 준선형 스케일링).

**★ 왜인가?**

```
포락선의 코어 변동: 진폭 33% / 섬 폭 ~15%
★ 그런데 검증(단일 공명)은 좁은 영역만 본다 → 포락선이 거의 상수 → 둔감

★ 진짜 결함은 경계에 있었다 (각평균 χ(ρ), m=12, κ⊥=1e-6):
  옛 자기장: 경계에서 χ ≈ 0.5~0.7   ← ★ T=0/1 경계로 평행 flux 유입
  새 자기장: 경계에서 χ ≈ 0.00      ← ★ Paul 포락선이 금지하는 것
```

> **★★ 검증이 통과했다고 모델이 맞는 것이 아니다.**
> **★★ 검증이 그 차이에 둔감했을 뿐일 수 있다.**

**★ 그리고 정합성 확인(A)의 반전:** 재구현한 converse-KAM/WBA **진단**은 원본과 일치했다
(검출 100%, 분류 93–96%). **결함은 진단이 아니라 자기장에 있었다.**

### ⑤ 도구 중복 / 데이터 소실 — A에서 잡음

**★ C3b 1차가 C2/C3a 원본을 fetch하지 않고 재구현했다.**
- 원본이 다른 브랜치에 있었고, **★ 원본도 격자 데이터를 .gitignore해 소실**시켰다.

**★ 대응 (C3b v2):**
- **★ 원본 코드를 fetch해서 썼다** (`consistency/orig/`, `bf5ca06` 핀 — `tools/`와 byte-identical).
- **★ 격자 데이터를 커밋한다** (`results/grid/*.npz`, .gitignore 예외).

> **★ 재현 불가능한 결과는 결과가 아니다.**

---

## 도구

| 파일 | 내용 | 검증 |
|---|---|---|
| `tools/field_kmm.py` | Kallinikos-MacKay 자기장 | 불변량 6×10⁻¹² 보존 |
| `tools/field_paul.py` | **★ Paul 자기장 (ψ(ψ−ψ̄) 포락선 포함)** | Chirikov S = 1.000, separatrix 1/8, 7/8 |
| `tools/converse_kam.py` | converse-KAM (MacKay 2018) | 섬 면적 0.3% 오차 |
| `tools/wba.py` | WBA (Duignan-Meiss 2023) | 게이트 통과 (섬 12–14 vs 바다 1) |
| `tools/residue.py` | Greene residue | R_O = +0.52, R_X = −0.69 |
| `vpd/` | V_PD (Paul 2022): `solver3d.py` + `paulsolve.py` | ε_crit=√κ⊥/2, √ε, 준선형 통과; smoke m=4 재현 |
| `results/grid/*.npz` | 격자 데이터 (t_c, dig, χ, T) | **★ 전부 커밋됨** |

**★ 소스 오브 레코드는 `tools/` 한 벌이다.** `consistency/orig/`는 그 byte-identical 핀이고,
`vpd/field.py::paul_field`·`vpd/diagnostics.py`는 SUPERSEDED 재구현본이다(사용 금지 — `consistency/SOURCE.txt`).
`tools/`의 도구는 QUASR/WBA 외부 원본이 아니라 논문에서 새로 구현했다(`tools/SOURCE.txt`).

---

## 문헌

- **MacKay (2018)** Reg. Chaotic Dyn. 23, 797 — converse-KAM 3D 이론
- **Kallinikos, MacKay & Martínez-del-Río (2023)** PPCF 65, 095021 — 토로이달 구현
  - 코드: https://github.com/dvmtz-1/cKAM
- **★ Paul, Hudson & Helander (2022)** JPP 88, 905880107 — V_PD, **그리고 이 프로젝트의 대상**
- **Duignan & Meiss (2023)** Physica D 449, 133749 — WBA
- **Vlad et al. (2002)** — λ_e, 자기장만으로 수송이 결정되지 않음

---

## ★ 아직 못 하는 것

**① cantori 가설은 여전히 가설이다.**
converse-KAM이 cantori를 못 보기 때문에 실패하는가?
관찰(m=36 카오스가 연결되어 있으나 전송에 10⁵ 회전 필요)은 **근거이지 증명이 아니다.**
**★ cantori의 turnstile flux를 직접 재야 한다.**

**② 섬 중심 foliation** (KMM §5). 섬의 elliptic 중심 foliation을 쓰면 converse-KAM이
섬과 카오스를 구분할 수 있다. **★ 그러면 위상 축이 완성된다.** 아직 안 했다.

**③ 유한 β.** 이 프로젝트는 전부 **주어진 자기장**에서의 분석이다.
**★ 압력이 자기장을 바꾸는 문제(유한 β)는 HINT이 필요하다.** 히로시마.

(전체 목록: `records/future_questions.md`)

---

## ★ 개발 이력 — 그리고 왜 이것이 중요한가

**★ C3b는 별도 브랜치에서 시작되었다. 그리고 그것이 사고를 낳았다.**

```
C3b 세션이 빈 트리에서 시작
  → ★ 원본 브랜치를 fetch하지 않음
  → ★ converse-KAM, WBA, 그리고 자기장을 재구현
  → ★ 자기장에서 Paul의 ψ(ψ−ψ̄) 포락선을 누락
  → ★ 그런데 검증을 전부 통과 (검증이 그 차이에 둔감)
  → ★ 정합성 확인(A)에서 자기장 불일치 발견 (진단은 멀쩡했다)
  → ★ 원본을 fetch해서 재실행 (C3b v2) — 결론은 유지·선명해짐
```

**★ 교훈 둘:**

> **① 도구는 원본이 하나다. 재구현하지 말고 fetch하라.**
> **② 격자 데이터를 커밋하라. 재현 불가능한 결과는 결과가 아니다.**

**★ 현재 상태 — 모든 작업이 이 브랜치(`claude/converse-kam-3d-y23ijj`)에 병합되어 있다.**
개발 이력은 `claude/v-pd-implementation-comparison-9aemtw` 브랜치에 보존.

---

## 재현

```bash
pip install numpy scipy pyamg

# C2: converse-KAM (원본)
python3 scripts/paul_ckam_campaign.py

# C3a: WBA (원본)
python3 scripts/wba_stage1.py

# C3b v2: V_PD + 세 축 비교 (Paul 명세 자기장, 원본 진단)
python3 -m vpd.validate1_paul     # Stage-1 게이트 (단일 섬)
python3 -m vpd.validate2_paul     # Stage-1 게이트 (카오스 층)
python3 -m vpd.stage2_paul        # 네 자기장: V_PD, ΔT
python3 -m vpd.stage3_paul        # 세 축 비교 + 상관
```

**★ 격자 데이터가 `results/grid/`에 커밋되어 있으므로, 재계산 없이 재분석 가능하다.**
(smoke test: `vpd.stage2_paul`의 m=4가 V_PD = 0.246 / 0.707 / 0.897을 이 환경에서 재현.)
