# ocean_box_model

[日本語](#日本語) | [English](#english)

---

# 日本語

## ocean_box_model

Toggweiler (1999) の海洋ボックスモデルを基礎として作成した、Fortran と Python による教育・研究用の海洋炭素循環モデルです。

本リポジトリでは、

- Fortran によるオリジナル実装
- Python による教育用実装
- Jupyter Notebook を用いた段階的な教材
- 演習問題・想定解答

を提供します。

Notebook は、海洋ボックスモデルを初めて学ぶ学生を対象に、

```
Energy balance model
    ↓
One-box
    ↓
Two-box
    ↓
Three-box
    ↓
Four-box
    ↓
Six-box
    ↓
Seven-box
```

という流れで、プログラミング・海洋循環・炭素循環を同時に学べるよう設計されています。

---

## Repository structure

```
ocean_box_model/

├── f90/
│   Fortran implementation
│
├── python/
│   Python implementation
│
├── notebook/
│   Bilingual (Japanese/English) teaching notebooks
│
├── notebook-solutions/
│   Solutions to notebook exercises
│
└── README.md
```

---

## Notebook contents

| Chapter | Contents |
|----------|----------|
| 00 | Energy balance model |
| 01 | One-box model |
| 02 | Two-box model |
| 03 | Three-box model |
| 04 | Four-box model |
| 05 | Six-box model |
| 06 | Seven-box model |

---

## Reference

Toggweiler, J. R. (1999), *Variation of atmospheric CO2 by ventilation of the ocean's deepest water*, Paleoceanography, **14**(5), 571–588. https://doi.org/10.1029/1999PA900033

---

# English

## ocean_box_model

**ocean_box_model** is an educational and research repository of ocean box models based on the six-box and seven-box models of Toggweiler (1999).

The repository contains

- original Fortran implementations,
- educational Python implementations,
- step-by-step Jupyter notebooks,
- exercise solutions.

The notebooks are designed for students with little background in either oceanography or scientific programming.

The learning sequence is

```
Energy balance model
      ↓
One-box Model
      ↓
Two-box Model
      ↓
Three-box Model
      ↓
Four-box Model
      ↓
Six-box Model
      ↓
Seven-box Model
```

Students gradually learn

- Python programming,
- matrix-based transport models,
- tracer transport,
- ocean circulation,
- carbon cycling,
- atmospheric CO2,
- ideal age,
- stable and radiocarbon isotopes.

---

## Repository structure

```
ocean_box_model/

├── f90/
│   Fortran source code
│
├── python/
│   Python source code
│
├── notebook/
│   Teaching notebooks (Japanese/English)
│
├── notebook-solutions/
│   Solutions to notebook exercises
│
└── README.md
```

---

## Reference

Toggweiler, J. R. (1999). *Variation of atmospheric CO2 by ventilation of the ocean's deepest water*. Paleoceanography, **14**(5), 571–588. https://doi.org/10.1029/1999PA900033
