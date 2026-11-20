# 01 — 배경지식 및 기존 실험·이미지 분석 맥락

> **문서 목적:** Claude가 이번 GFAP fluorescence microscopy 데이터 분석의 생물학적·실험적 배경을 이해하도록 하는 참고 문서. **현재 코딩 과제의 직접 명세는 `03_Coding_Task_and_Analysis_Spec.md`를 우선**한다.
>
> **기준 시점:** 2026-11-19. 이 문서에는 사용자가 제공한 Rose 이메일 원문, 실험 사실, 대화에서 도출한 추론을 혼합하지 않도록 출처 성격을 구분했다.
>
> **증거 구분:** `[ROSE 원문]` 직접 인용 또는 제공된 이메일 화면에서 읽어 옮긴 내용; `[사용자/미팅 정리]` 사용자가 제공한 기록·사진·미팅 요약; `[해석/권고]` 후속 작업을 위해 제안하는 분석 방법이며 Rose가 직접 지시한 것으로 간주하면 안 됨.

## 1. 연구 환경과 핵심 개념

- **연구실:** McGill University, Ribeiro-da-Silva (RDS) Laboratory. 지도교수 Prof. Alfredo Ribeiro-da-Silva, 연구 멘토 Rose Rodrigues.
- **현재 연구:** PHAR 599 D1/D2에서 mouse spinal cord 조직의 **GFAP immunofluorescence staining**, 항체 비교, 두께 방향(signal penetration) 문제의 최적화, 이후 astrocyte image annotation/segmentation에 사용할 데이터 확보.
- **주요 조직:** free-floating 척수 횡단 절편. 실험 프로토콜상 **50 µm 절편**이 자주 사용됨. *명목상 절편 두께 50 µm가 모든 `.lif`의 실제 촬영된 z-depth와 같다는 뜻은 아니다.* 실제 z-depth는 원본 메타데이터로 확인한다.
- **GFAP (glial fibrillary acidic protein):** astrocyte intermediate filament 관련 표지. 세포 돌기(process)를 보여주지만 GFAP-positive 구조만으로 개별 astrocyte의 전체 세포막이나 각 돌기의 소속을 완벽히 구별할 수 있는 것은 아님.
- **Primary antibody:** Mouse anti-GFAP / Rabbit anti-GFAP / Chicken anti-GFAP. host species와 항원 target을 구별할 것.
- **Secondary antibody:** primary host에 결합하며 형광 표지가 달림. 예: Donkey anti-Mouse Alexa Fluor 647(AF647). 현재 농도곡선 실험은 **새 Mouse anti-GFAP primary**와 **anti-Mouse AF647 secondary**의 조합.
- **DAPI:** DNA/핵 표지. GFAP signal의 정량을 DAPI 채널과 혼동하지 말 것.
- **z-stack:** 같은 시야를 조직 깊이(z) 여러 위치에서 촬영한 3D 이미지. 통상 한 GFAP 채널의 shape는 `(Z, Y, X)`로 정리할 수 있지만, `.lif`의 실제 축 순서는 읽기 라이브러리/메타데이터에 따라 다름. **확인 없이 reshape하지 말 것.**
- **Maximum-intensity projection (MIP, max projection):** 동일 (x,y)에서 z축 최대 intensity를 취한 이미지. 높은 intensity를 보이는 process 관찰에는 유용하나, 평균 밝기는 **z-slice 개수**나 depth coverage 차이의 영향을 받을 수 있음.
- **Pixel bit value:** 이미지 픽셀에 저장된 디지털 intensity 값. 16-bit 원본이라면 일반적으로 0~65535 범위의 정수일 수 있으나 센서/채널 획득 설정을 확인해야 함. **bit depth와 실제 측정 단위(photons, counts, intensity)를 혼동하지 말 것.**
- **Autofluorescence / background:** GFAP-specific staining과 무관한 신호. 조직 자체 배경, perfusion 문제, 비특이 결합, acquisition noise 등이 가능한 원인. 원인을 단정해서는 안 됨.
- **Deconvolution:** optics의 PSF 영향을 계산적으로 줄여 구조를 더 선명하게 보이도록 하는 영상 처리. 이번 파일에서는 이름에 `_lng`가 붙은 것이 **deconvolved + full-16-bit-range normalized** 버전이라고 Rose가 설명했다. 원본과 같은 절대 intensity 척도로 비교하면 안 됨.

## 2. 왜 깊이에 따른 GFAP signal이 중요한가?

### 관찰된 문제

과거 Mouse/Rabbit/Chicken 각각을 비교한 staining에서 모든 antibody가 형광 signal을 냈지만, **stack 중심부 signal이 표면보다 약한 현상**이 있었다. 일부 antibody에서 특히 심했다. 기술적으로는 다음이 모두 가능함.

1. Primary나 secondary antibody가 두꺼운 조직 내부까지 불균일하게 침투함.
2. 깊이에 따른 광학적 산란, 흡수, imaging/detection 차이.
3. 조직 상태와 배경/자체형광 차이.
4. 깊이에 따른 GFAP 구조 분포의 실제 생물학적 차이.

**중요:** z-middle이 어둡다는 관찰은 *penetration과 양립하지만*, 이미지만으로 primary/secondary 침투 문제를 확정하는 증거는 아님. 그래서 incubation time, fluorophore 등의 변수를 통제하는 troubleshooting이 진행되었다.

### 왜 AF647로 통일했나?

Rose가 명시한 의도는 (i) far-red detector의 spectral acquisition window 이점, (ii) 더 긴 파장이 두꺼운 조직에서 일반적으로 산란이 적을 수 있는 점, (iii) Mouse/Rabbit/Chicken primary 조건을 보다 공정하게 비교할 수 있는 점이다. **다른 fluorophore에서 나온 수치끼리 raw intensity를 비교하면 fluorophore/detector confounding이 발생**한다. AF647을 맞춰도 antibody affinity, epitope, secondary labelling efficiency가 동일해지는 것은 아님.

## 3. 이전에 Rose가 제공한 두 PDF의 역할

두 PDF는 **과거 Mouse/Rabbit/Chicken single-antibody comparison의 결과물**이다. **현재 1:250–1:2000 Mouse-only concentration curve의 결과물이라고 착각하지 말 것.** Claude에게 두 PDF를 추가로 주면 분석의 예시 시각자료로 사용할 수 있다.

### `GFAP_maxproj_all_magma.pdf` (3 pages)

- p.1 Mouse anti-GFAP triplicates (`Ms_1`, `Ms_2`, `Ms_3`): 63×, raw GFAP 638 nm, z-max projection, magma display, per-image 99.8th-percentile 표시 범위, pixel-intensity histogram.
- p.2 Rabbit triplicates (`Rb_1`, `Rb_2`, `Rb_3`), 같은 유형.
- p.3 Chicken triplicates (`Ck_1`, `Ck_2`, `Ck_3`), 같은 유형.
- 보여주는 것: GFAP morphology의 전체 투영, 상대적으로 강/약한 staining, histogram shape.
- **주의:** 이미지마다 표시의 상한값이 다르게 설정된 경우, 눈으로 보이는 brightness 차이는 raw intensity 차이와 1:1 대응하지 않음. raw values/table도 병기할 것.

### `GFAP_depth_slices_magma.pdf` (2 pages)

- p.1 Mouse/Rabbit/Chicken 각각의 대표 stack에서 0%, 25%, 50%, 75%, 90% depth slice 예시와 slice별 99.8th-percentile photon-count depth curve.
- p.2 각 antibody triplicates의 **per-slice 99.8th-percentile GFAP intensity vs depth** 선 그래프.
- 보여주는 것: z-depth에 따른 signal loss 또는 유지, profile shape, replicate간 variability.
- p.1 표기된 깊이는 시료별 약 20~25 µm 수준까지 표시된다. **이것을 곧바로 section 두께 50 µm와 일치한다고 가정하지 말 것.**

### 99.8th percentile에 관한 해석 주의

Rose는 아주 밝은 소수 픽셀(outlier)이 표시 상한을 지배하는 문제를 피하려고 **99.8th percentile을 `normalization max`**로 사용했다고 설명했다.

- **표시/시각화:** 예를 들어 한 이미지의 99.8th percentile을 display vmax로 정하는 것은 타당한 예시.
- **정량:** raw 평균 intensity를 계산하려는 목적이라면 **먼저 raw 데이터에서 mean을 계산**하고, percentile-normalized display값과 별도로 보고한다.
- 두 operation을 혼동해서 `mean(raw / each_image_p99.8)`만 보고 서로 다른 농도의 절대 signal을 비교하면, 원래 농도에 따른 밝기 차이를 지울 위험이 있다.
- `p99.8`, raw mean, background-corrected raw mean, display vmax를 **서로 다른 column**으로 기록할 것.

## 4. 과거 사용자의 이미지 분석 경험 (참고; 이번 과제와 혼동 금지)

[사용자/기존 연구 기록] 사용자는 RDS 연구실의 과거 프로젝트 및 별도 이미지 분석에서 Fiji/ImageJ를 사용해 ROI masking, 채널별 처리, stack indexing, histogram 생성, CSV export, ratiometric fluorescence 정량화 등을 다뤘다. 예전 코드의 구현 언어는 주로 Fiji/ImageJ macro였을 수 있고, 지금 Rose는 **Python으로 `.lif`를 읽고 정량 분석을 처음부터 독립적으로 구현**해보라는 취지다.

별도의 예전 astrocyte workflow에서는 tdTomato/GFAP/Sox9, threshold-based pre-mask, napari correction을 다뤘다. 이 과거 경험은 이번 코드 개발에 도움이 되지만, **roGFP2 ratio 분석·astrocyte segmentation·이번 Mouse-GFAP concentration curve는 서로 동일한 데이터셋/측정값이 아니다.**

## 5. Rose 이메일 원문 ① — GFAP AF647 penetration troubleshooting (사용자가 채팅에 직접 제공한 영어 본문)

> 아래는 사용자가 직접 붙여넣은 Rose 이메일을 보존한 것. 공백·강조 마크다운은 전달 편의를 위해 정리했지만 영문 내용은 원문을 유지했다.

```text
Hi Sewon,

I've been imaging the samples, and while all the antibodies have signal, we're having some penetration issues (low signal in the middle). Let's start off by testing out 3 conditions, all stained at the same time:

1. normal secondary antibody protocol (2hrs room temp)
2. secondary antibody for 4hrs room temp
3. secondary antibody overnight in the cold room (this stain would finish a day after the others, so the sections should go on a different plate).

For these, let's do just the single antibody stains (Ms, Rb, Ch, separately). And this time, to control for as many 'factors' as possible, all the stains get AlexaFluor 647. I believe we have all these secondaries in Donkey, but double check before starting the stains in case you need to use Goat for blocking and primary antibody solutions. PS: make sure to include no primary controls, but we really only need like, 2 sections of spinal cord in any of these controls.

Reasoning 1: when running multiple excitation wavelengths in the same imaging track, the only detector that is guaranteed to not be cut off early is the far red, and since we always want the most detail in GFAP, that's why I always put GFAP in AF647 in 'real' experiments.

Reasoning 2: in microscopy, longer wavelengths give better imaging in 'thick' tissue because the longer the wavelength, the less scatter. This logic is most obvious when you employ a technique like two-photon imaging, which employs full infrared. but nonetheless, we have a better chance to minimize low-signal-in-the-middle by taking advantage of the longest possible wavelength we have access to.

Reasoning 3: Choosing the same excitation / emission parameters, we can more directly compare photon counts between the 3 antibodies. Yes, technically they are "different" antibodies since the antigen is different, but it's the closest we can get.

As it stands, Ch-GFAP is looking the best, even though it has the highest background signal. But until we do these tests, I can't say for sure if it's because the antibody is superior or because it's the one that has AF647.

Enjoy your break!
```

## 6. Rose 이메일 원문 ② — 과거 antibody별 depth-profile 분석 결과 (사용자가 채팅에 직접 제공)

```text
Hi Sewon,

Forgot to book the microscope this week, but at least I can show you the data I have from the first stain.

Triplicates of each individual antibody stain, and Rabbit looks much worse than the other two on max projection. So then I computed the change in intensity across the depth of the stack for each (I use 99.8th percentile as a normalization value, so “not too many” pixels are saturated but still none of the data is skewed by a random much-brighter-than-normal pixel), and clearly there’s just about zero signal in the Rabbit in the middle section, while mouse and chicken are more stable.

Most antibodies tend to have “more” intensity in the top quarter of the stacks, this is usually the case in thick tissue, so I’m not too worried about how both mouse and chicken have a higher hump on the top 20% of the stack.

We can do actual numbers and computing later on, but this already gives you an idea of the first experimental results you have.

I’ll let you know when I finish the imaging next week!
```

## 7. Claude에게 넘길 때 체크할 사항

- 먼저 **현재 Mouse-only concentration curve**인지, **과거 3-antibody 비교**인지 프로젝트/데이터를 구분한다.
- `GFAP_depth_slices_magma.pdf`, `GFAP_maxproj_all_magma.pdf`는 분석 **예시 자료**이지 새로운 `.lif`에 대한 결과가 아니다.
- Rose 원문에 명시된 주장과 이 문서의 통계/코딩 제안을 혼동하지 않는다.
- 데이터나 acquisition metadata 없이 실제 intensity values, z-depth, signal/background 수치를 만들어내지 않는다.
- 이후 코드 명세는 반드시 `03_Coding_Task_and_Analysis_Spec.md`를 우선 참고한다.
