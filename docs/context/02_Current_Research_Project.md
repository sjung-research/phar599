# 02 — 현재 프로젝트의 전체 목적, 연구 흐름, 실제 진행 상황

> **Claude를 위한 목적:** 지금 맡길 `.lif` intensity 분석이 *어떤 연구 질문*을 해결하고, 최종 astrocyte segmentation 연구에 어떻게 이어지는지 이해하기 위한 프로젝트 지도.
>
> **상태 기준:** 2026-11-19. **이번 주에 수행해야 하는 코딩만** 필요하면 `03_Coding_Task_and_Analysis_Spec.md`를 먼저 읽어도 된다.
>
> **출처 라벨:** `[ROSE 원문]` 사용자가 제공한 이메일/스크린샷; `[미팅 요약]` 사용자가 제공한 이전 Rose 미팅의 2차 요약; `[계획/추론]` 연구 목적을 일관된 구조로 묶기 위한 제안. **추론을 Rose가 확정한 일정·방법으로 바꾸어 읽지 말 것.**

## 1. 한눈에 보는 연구 목표

### Ultimate objective (최종 지향점)

**척수 조직에서 개별 protoplasmic astrocyte의 복잡한 morphology를 일관되게 시각화·annotation·segmentation하여, 나중에 형태학적 정량분석이 가능한 reproducible imaging-to-analysis pipeline을 구축한다.**

구조:

```text
Mouse spinal cord specimens
   ↓
GFAP antibody / staining / imaging optimization
   ↓
Multimodal confocal imaging (GFAP + sparse tdTomato + Sox9)
   ↓
Ground-truth labels: initial masks + manual correction (napari)
   ↓
Curated dataset (existing ~50 images + future images)
   ↓
Train/adapt/evaluate computer-vision segmentation
   ↓
Quantify single-astrocyte morphology
```

**주의:** 이 전체 파이프라인은 연구실의 장기 목표를 정리한 것이다. 현재 학부 프로젝트에서 모든 단계를 이미 완료했다거나 모든 후속 기술적 선택이 확정됐다는 뜻은 아니다. 현재 실제 수행한 작업의 상당 부분은 staining troubleshooting과 concentration curve이다.

## 2. 연구실 및 사용자 역할

- 학생: Sewon (McGill PHAR 599 D1/D2); 지도교수 Prof. Alfredo Ribeiro-da-Silva, 직접 멘토 Rose Rodrigues.
- Rose가 새 antibody order, microscopy imaging, 일부 Python 분석과 annotation workflow(기존 napari 환경)를 주도하고, 사용자가 staining 단계 및 일부 quantitative image-analysis task를 맡아왔다.
- 새 Mouse anti-GFAP primary의 `1:250, 1:500, 1:1000, 1:2000` **serial dilution을 실제로 만드는 단계는 Rose가 tube 간 transfer/vortex/mixing을 직접 수행했다**고 사용자가 설명했다. 학생은 Day 2 Donkey anti-Mouse AF647 secondary/DAPI/wash를 스스로 완료했다. 따라서 이번 intensity 결과의 불확실성을 학생의 primary pipetting error로 단정해서는 안 된다.
- Rose가 먼저 사용자가 독립적으로 Python 코드를 설계한 뒤 자신의 이전 분석 Python script와 **방법을 비교**하고 싶어 한다. 정답 코드가 반드시 이미 공유된 것은 아님.

## 3. Rose 미팅의 원래 프로젝트 목표: 왜 좋은 GFAP staining이 필요한가

[미팅 요약] 2026년 9월 말 Rose와의 약 70분 미팅을 사용자가 제공한 요약에 따르면:

- **Sparse tdTomato astrocytes + GFAP** 이미지로 protoplasmic astrocyte segmentation model의 training data를 만드는 것이 원래 프로젝트의 장기 목표다.
- **Sox9**는 astrocyte nucleus 또는 cell centre를 식별할 수 있는 보조 표지다.
- 이전에 약 **50개의 tdTomato + GFAP 이미지**를 취득했다.
- 전부 수동으로 astrocyte 구조를 delineate하려면 매우 오랜 시간이 걸렸다(미팅 요약에 따르면 세포 하나당 약 3일 사례).
- Rose는 thresholding, 특히 **triangle-type threshold 기반 pre-mask**를 만들고 이를 사람이 수동 교정하는 방식을 시도했다. **napari**에서 기존 label을 수정하며 맞춤형 hotkey가 있다.
- 예전 Mouse anti-GFAP 제품/lot 문제 때문에 기존 50장 labeling이 약 3–4장 정도 진행된 시점에서 중단되었다.
- 대안 Mouse / Rabbit / Chicken anti-GFAP를 검사했지만 가장 밝은 영상이 반드시 segmentation에 가장 유리한 것은 아니다. **background, 복잡한 GFAP process overlap, 세포 간 경계 해석 가능성**이 중요하다.
- Mouse primary 중 기존 데이터와 일치하는 제품은 비쌌지만 original 50 images와 staining 특성을 맞출 수 있어 continuity가 크다.
- napari를 통해 기존 이미지 annotation이 재개되면, 이후 추가 조직 획득/이미징과 model-development로 넘어갈 수 있다. Cellpose 관련 협업 가능성은 당시 논의된 장기 구상이지 확정된 공동연구나 선정된 모델이 아니다.

### 연구 Aim을 4개로 묶으면

| Aim | 핵심 | 프로젝트 단계/현황 |
|---|---|---|
| 1 | GFAP antibody, AF647 imaging, secondary incubation/penetration, Mouse concentration optimization | **진행 중**; 일부 비교 이미지 획득; 최종 concentration 미확정 |
| 2 | tdTomato/GFAP/Sox9 multimodal spinal-cord image dataset | 기존 ~50장 있음; 새 샘플 획득은 추후 |
| 3 | threshold pre-mask + napari manual correction으로 ground-truth annotation | 중단된 기존 작업 재개 계획; 환경 설치 문제 발생 |
| 4 | Automated single-astrocyte segmentation 및 morphology quantification | broader downstream goal; **이번 `.lif` 농도코드의 직접 목표는 아님** |

## 4. 지금의 집중 과제: Mouse anti-GFAP concentration curve

### 왜 시행했나?

기존에 좋은 morphology를 보였던 Mouse GFAP 항체와 일치하는 제품이 새로 도착했다. **새 aliquot의 최적 working concentration**을 알아보기 위해 primary antibody만 2-fold로 serial dilute했다.

- Primary Mouse anti-GFAP final dilution: `1:250 → 1:500 → 1:1000 → 1:2000`.
- 고정: Donkey anti-Mouse AF647 secondary; DAPI; 동일 실험 protocol에서 가능한 한 일관된 acquisition.
- 현재 `.lif` 원본에는 **각 농도당 이미지 하나와 그에 대응하는 `_lng` deconvolved image 한 버전**이 있다고 Rose가 명시했다.
- **현재 우선 질문:** four concentrations의 signal/contrast/background/depth retention을 정량화해 최적 농도를 판단할 수 있는가?
- Rose의 2026-11-13 결과: 네 aliquot의 사진은 4개 농도에서 거의 같아 보였고, Emma의 샘플은 농도별 intensity 차이를 보였다. **staining itself worked**, 그러나 **final optimal concentration은 아직 결정하지 않음**.
- Rose는 ipsilateral (`ipsi`) / contralateral (`contra`) 촬영도 추가하여 실제 정량화하려 했다. 현재 SharePoint로 어떤 추가 이미지가 전달되는지는 파일을 열어 확인할 것. 이번 `.lif`에 양쪽 정보가 모두 있다고 미리 가정하지 않는다.

### 절대 혼동하면 안 되는 것

1. **이번 concentration curve**는 primary **Mouse 한 종류의 희석비 비교**이다.
2. **과거 3-antibody comparison**은 Mouse/Rabbit/Chicken primary를 서로 비교했다. 관련 예시 PDF 두 개의 내용은 과거 비교다.
3. **과거 50장 tdTomato/GFAP dataset**은 napari annotation 목적의 이미징 데이터다. 이번 원본 `.lif` 농도시험과 동일 데이터라는 근거는 없다.
4. `_lng` version은 intensity comparison용이 **아니다**. Rose가 전체 16-bit range normalization을 명시했기 때문이다.

## 5. 2026-11-19 기준 프로젝트 진행 현황

| 작업 | 상태 | 다음 단계 |
|---|---|---|
| 새 Mouse GFAP 확보 | 완료 | 농도 선택 검증 |
| Mouse 2-fold concentration staining | 완료, signal 있음 | raw `.lif` 정량 분석 |
| 농도별 intensity에 따른 최종 concentration | 미확정 | raw MIP/depth/background/histogram 분석 + Rose 확인 |
| Original 50-image annotation | 시작/재개 희망 | napari environment 확보 |
| Lab computer의 Rose napari/conda environment 복제 | Rose가 시도했지만 conda 설치 문제 | 이후 재시도 필요 |
| 이번 주 학생의 Python 분석 | 신규 요청됨 | SharePoint에서 `.lif` 접근 후 독립 구현 |
| 다음 specimen 획득 | 준비 단계 | Rose tamoxifen injections, tissue collection 예정(몇 주 이상 지연 가능) |

Rose는 11월 중순 집 이사 준비와 Tesla의 추가 수의사 방문 일정이 있어 당장 lab 환경 구성을 지속하기 어려웠다. **이것은 배경 일정일 뿐 코딩에 필요한 분석 변수는 아니다.**

## 6. Rose 이메일 원문 ① — 새 Mouse antibody 도착 / 농도곡선 제안 (10월 30일경 스크린샷)

> 화면에 보인 영문 본문을 가급적 그대로 전사. 스타일 및 날짜 표기는 별도 메타데이터이며, 원문에서 수의사/이사 이야기는 없음.

```text
Hi Sewon,

Alfredo cleared the order for the Mouse GFAP that matches my previous data, and the vials have arrived!

Next week, you could run a concentration curve of “your” aliquot. It’ll be 1:250 to 1:2,000 in a serial dilution, but I don’t remember if you’ve done one of those before? Let me know if you have that protocol already, in which case you’re good to do the stain (AF647 + DAPI on day 2) all on your own! Or if you’d like to check with me during the day 1, I’ll be in lab every day but Thursday, from 3-6pm.

Let me know!
```

## 7. Rose 이메일 원문 ② — 농도곡선 결과 및 napari 재개 계획 (2026-11-13, 제목: `Project updates!`)

```text
Hi Sewon,

I took a look at your slides, and somehow, your aliquot of mouse anti GFAP looks exactly the same in all concentrations (unlike Emma’s, which changes intensity based on concentration). I want to go back and take more images from both ipsi and contra sides, to get actual quantification. So I don’t have a “final concentration” for you yet, though the stain absolutely worked.

BUT! We’re just about ready to start staining new samples. I just need to 1) do tamoxifen injections and 2) collect tissue. This will take a couple of weeks minimum, probably more because I’m currently in the middle of packing up to move apartments and poor Tesla still has more upcoming vet appointments.

I think in the meantime, since we’re keeping the same primary antibody as the original experiment, we can actually get into labelling the data from my original 50 images!!!!!!!!!!

On my calendar, I put in next Tuesday to (attempt to) set up an exact replica of my Napari environment on the lab computer, so that I can teach you how to open the image files and use the hotkeys I created in order to fix the labels on the images. If I’m successful with that Tuesday, I could do Wednesday or Friday to show you how it all works. How does that sound?

Sorry for the long email!
```

## 8. Rose 이메일 원문 ③ — conda 문제로 분석을 먼저 맡김 (2026-11-18, 제목: `I failed :o/`)

```text
I'm losing the fight with conda, and I need to finish packing my house for my move Saturday so I have to tap out for this week. In the interim since I think Friday is no-go, how about I give you access to the .lif file for the images that I took of the concentration curve, and you try creating/running the analysis of the bit values and variance over stack yourself?

Are you okay with playing around with that this week while I'm moving houses and getting set up? We can let you try it on your own first, then we can compare python files of how I did it last time vs what you make on your own!

Let me know what you think!
```

## 9. Claude에게 꼭 전달할 결론

- **이 코딩은 완성된 GFAP segmentation model을 당장 만들라는 과제가 아니다.** 이미지의 raw fluorescence signal을 분석하여 농도 선택 문제에 답하는 독립적 정량화 작업이다.
- 분석 결과가 충분하지 않으면 **"optimal concentration cannot yet be determined"**라고 솔직히 결론 내릴 수 있어야 한다.
- Rose의 예전 Python code는 사용자가 독립적인 해결안을 만든 후 비교하는 대상으로 남겨두는 것이 의도에 맞다.
- 이미지를 볼 때 예전 GFAP 3-antibody PDF를 **시각화/통계 설계 참고**로만 쓰고, 새 실험의 수치를 거기서 복사해서는 안 된다.
- 다음 문서 `03_Coding_Task_and_Analysis_Spec.md`가 구현에 관한 우선 명세다.
