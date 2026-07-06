# 00 モデルとプログラミング（EBM）

## 課題 1 / Exercise 1

### 問題 / Question

λ を変えると ECS はどう変わるか。  
How does ECS change when λ is varied?

### 想定解答 / Expected Answer

ECS は

\[
ECS=\frac{F_{2\times CO_2}}{\lambda}
\]

で決まる。

したがって、λ が小さいほど平衡気温上昇は大きくなる。

つまり、気候フィードバックが弱いモデルほど温暖化しやすい。

ECS is determined by

\[
ECS=\frac{F_{2\times CO_2}}{\lambda}.
\]

Therefore, a smaller λ leads to a larger equilibrium temperature increase. In other words, a model with weaker climate feedback exhibits greater warming.

---

## 課題 2 / Exercise 2

### 問題 / Question

熱容量 C を変えると何が変わるか。  
What changes when the heat capacity C is varied?

### 想定解答 / Expected Answer

熱容量が大きいほど温度変化はゆっくりになる。

一方で平衡状態は

\[
T=\frac{F}{\lambda}
\]

なので、最終温度は変化しない。

つまり、C は応答速度を決めるパラメータである。

A larger heat capacity slows the temperature response.

However, the equilibrium temperature is

\[
T=\frac{F}{\lambda},
\]

so the final equilibrium temperature does not change.

Thus, C controls the response timescale rather than the equilibrium state.

---

## 課題 3 / Exercise 3

### 問題 / Question

なぜ TCR は ECS より小さいか。  
Why is TCR smaller than ECS?

### 想定解答 / Expected Answer

CO2 が増えている途中では海洋が熱を吸収し続けるため、まだ平衡状態に達していない。

そのため、約70年後の温暖化（TCR）は、十分長時間後の温暖化（ECS）より小さい。

During increasing CO2, the ocean continues to absorb heat, so the climate system has not yet reached equilibrium.

Therefore, the warming after about 70 years (TCR) is smaller than the long-term equilibrium warming (ECS).

---

## 課題 4 / Exercise 4

### 問題 / Question

Mini CMIP で ECS と TCR が一致しない理由を説明せよ。  
Explain why ECS and TCR are different in the Mini CMIP experiment.

### 想定解答 / Expected Answer

ECS は主に気候フィードバック λ によって決まる。

一方、TCR は λ に加えて海洋熱容量 C の影響も受ける。

そのため、ECS が同程度でも、海洋の熱吸収が大きいモデルでは TCR は小さくなる。

ECS is mainly determined by the climate feedback parameter λ.

In contrast, TCR is influenced by both λ and the ocean heat capacity C.

Therefore, models with similar ECS values can have different TCR values depending on the strength of ocean heat uptake.

---

## 課題 5 / Exercise 5

### 問題 / Question

EBM と海洋ボックスモデルに共通するプログラム構造を説明せよ。  
Explain the common program structure shared by the EBM and ocean box models.

### 想定解答 / Expected Answer

どちらも

1. 現在の状態を用いて
2. 各変数の変化量を計算し
3. 次の時刻の値を計算し
4. 状態を更新する

という時間積分を繰り返している。

扱う変数は異なるが、プログラムの基本構造は同じである。

Both models repeatedly perform time integration by

1. using the current state,
2. calculating the tendency of each variable,
3. computing the values at the next time step, and
4. updating the model state.

Although the variables differ, the overall program structure is the same.

---

# 01 One-box model

## 課題 1 / Exercise 1

### 問題 / Question

Export Production（CEPS）を増やすと PO4、DIC、O2 はどう変化するか。  
How do PO4, DIC, and O2 change when the export production (CEPS) is increased?

### 想定解答 / Expected Answer

輸出生産が強くなるため、

- PO4 は減少する。
- DIC は減少する。
- O2 は増加する。

また、それぞれの変化速度も大きくなる。

As export production becomes stronger,

- PO4 decreases.
- DIC decreases.
- O2 increases.

The rates of change also become larger.

---

## 課題 2 / Exercise 2

### 問題 / Question

One-box モデルの最大の欠点は何か。  
What is the main limitation of the one-box model?

### 想定解答 / Expected Answer

生物ポンプによって表層から取り除かれた物質の行き先を表現できないことである。

実際には沈降粒子は深層へ輸送されるが、One-box モデルではその過程を表現できない。

The main limitation is that the model cannot represent the fate of materials removed from the surface by the biological pump.

In reality, sinking particles are transported to the deep ocean, but this process cannot be represented in the one-box model.

---

## 課題 3 / Exercise 3

### 問題 / Question

なぜ Two-box モデルが必要なのか。  
Why is a two-box model necessary?

### 想定解答 / Expected Answer

現実の海洋では、

- 表層で栄養塩が消費される。
- 深層で有機物が分解され、栄養塩が再無機化される。

という循環が存在する。

この上下方向の輸送を表現するためには、少なくとも表層と深層の2つのボックスが必要となる。

In the real ocean,

- nutrients are consumed in the surface layer, and
- organic matter is decomposed in the deep ocean, where nutrients are remineralized.

To represent this vertical transport and recycling, at least two boxes (surface and deep ocean) are required.

---
# 02 Two-box model

## 課題 1 / Exercise 1

### 問題 / Question

表層と深層で PO4 の変化が逆になる理由を説明せよ。

Explain why PO4 changes in opposite directions in the surface and deep boxes.

### 想定解答 / Expected answers

表層では植物プランクトンが PO4 を消費する。

沈降粒子は深層で分解されるため、深層では PO4 が再無機化される。

そのため、

- 表層では PO4 が減少する。
- 深層では PO4 が増加する。

In the surface box, phytoplankton consume PO4.

Sinking organic particles are decomposed in the deep ocean, where PO4 is remineralized.

Therefore,

- PO4 decreases in the surface box.
- PO4 increases in the deep box.

---

## 課題 2 / Exercise 2

### 問題 / Question

循環速度 T を大きくすると何が起こるか。

What happens when the circulation rate T is increased?

### 想定解答 / Expected answers

表層と深層の交換が速くなる。

その結果、

- 表層と深層の濃度差は小さくなる。
- 系全体がより均一な状態に近づく。

Increasing T accelerates the exchange between the surface and deep boxes.

As a result,

- the concentration difference between the two boxes becomes smaller, and
- the system approaches a more homogeneous state.

---

## 課題 3 / Exercise 3

### 問題 / Question

Export Production を増やすと DIC、PO4、O2 はどのように変化するか。

How do DIC, PO4, and O2 change when export production is increased?

### 想定解答 / Expected answers

表層では生物生産が強くなるため、

- PO4 は減少する。
- DIC は減少する。
- O2 は増加する。

一方、深層では有機物が分解されるため、

- PO4 は増加する。
- DIC は増加する。
- O2 は減少する。

As biological production becomes stronger in the surface ocean,

- PO4 decreases.
- DIC decreases.
- O2 increases.

In the deep ocean, organic matter is decomposed, so

- PO4 increases.
- DIC increases.
- O2 decreases.

---

## 課題 4 / Exercise 4

### 問題 / Question

One-box と Two-box の違いを説明せよ。

Explain the difference between the one-box and two-box models.

### 想定解答 / Expected answers

One-box モデルは保存則を学ぶための最も単純な教育モデルであり、物質輸送を表現できない。

Two-box モデルでは表層と深層を分けることで、

- 生物ポンプ
- 深層での再無機化
- 鉛直循環

を表現できるようになる。

The one-box model is the simplest educational model for learning conservation laws and cannot represent material transport.

By separating the surface and deep oceans, the two-box model can represent

- the biological pump,
- remineralization in the deep ocean, and
- vertical circulation.

---

# 03 Three-box model

## 課題 1 / Exercise 1

### 問題 / Question

高緯度表層と低緯度表層を分ける理由を説明せよ。

Explain why the surface ocean is divided into high-latitude and low-latitude boxes.

### 想定解答 / Expected answers

現実の海洋では、

- 深層水形成
- 栄養塩濃度
- CO2 の海洋・大気交換

が高緯度と低緯度で大きく異なる。

その違いを表現するため、高緯度表層と低緯度表層を別々のボックスとして扱う。

In the real ocean,

- deep-water formation,
- nutrient concentrations, and
- air–sea CO2 exchange

differ substantially between high and low latitudes.

To represent these differences, the surface ocean is divided into separate high-latitude and low-latitude boxes.

---

## 課題 2 / Exercise 2

### 問題 / Question

High-latitude box が大気 CO2 に重要な理由を説明せよ。

Explain why the high-latitude box is important for atmospheric CO2.

### 想定解答 / Expected answers

高緯度では深層水形成が起こる。

そのため、

- 深層へ炭素を輸送する。
- 深層から CO2 を大気へ放出する。

という両方の役割を担っている。

このため、高緯度ボックスの状態は大気 CO2 濃度に強く影響する。

Deep-water formation occurs at high latitudes.

Therefore, the high-latitude box plays both of the following roles:

- transporting carbon into the deep ocean, and
- releasing CO2 from the deep ocean to the atmosphere.

As a result, the state of the high-latitude box strongly influences atmospheric CO2 concentration.

---

## 課題 3 / Exercise 3

### 問題 / Question

Two-box では再現できず、Three-box で再現できることは何か。

What can be represented in the three-box model but not in the two-box model?

### 想定解答 / Expected answers

高緯度と低緯度の役割分担を表現できることである。

例えば、

- 高緯度での深層水形成
- 高緯度でのガス交換
- 低緯度での生物生産

を別々に扱うことができる。

The three-box model represents the different roles of high and low latitudes.

For example, it separately represents

- deep-water formation at high latitudes,
- air–sea gas exchange at high latitudes, and
- biological production at low latitudes.

---

## 課題 4 / Exercise 4

### 問題 / Question

Three-box モデルの限界は何か。

What is the main limitation of the three-box model?

### 想定解答 / Expected answers

高緯度海洋を一つのボックスで表しているため、

- 北大西洋
- 南大洋

を区別できない。

実際には両者は異なる循環構造を持つため、この違いは表現できない。

The high-latitude ocean is represented by a single box, so the model cannot distinguish between

- the North Atlantic and
- the Southern Ocean.

In reality, these regions have different circulation systems, and this difference cannot be represented.

---

## 課題 5 / Exercise 5

### 問題 / Question

Four-box モデルでは何を改善したいと思うか。

What improvements would you make in the four-box model?

### 想定解答 / Expected answers

Three-box モデルの高緯度ボックスを

- 南大洋
- 北大西洋

に分離することで、

深層水形成と深層水の湧昇を別々に表現できるようになる。

これにより、現実の全球熱塩循環や Toggweiler (1999) の考え方により近いモデルとなる。

The high-latitude box in the three-box model is separated into

- the Southern Ocean and
- the North Atlantic,

allowing deep-water formation and deep-water upwelling to be represented separately.

This makes the model more consistent with the real global thermohaline circulation and the conceptual framework of Toggweiler (1999).

---

# 04-01 Four-box model

## 課題 1 / Exercise 1

### 問題 / Question

3-box モデルで高緯度を 1 つの箱にした場合、何が表現できないか。  
What cannot be represented when high latitudes are treated as one box in the 3-box model?

### 想定解答 / Expected answers

3-box では、高緯度が 1 つの箱なので、北大西洋の沈み込みと南大洋の湧昇を区別できない。  

In the 3-box model, high latitudes are represented by one box, so North Atlantic sinking and Southern Ocean upwelling cannot be distinguished.

---

## 課題 2 / Exercise 2

### 問題 / Question

4-box モデルで N box を追加する科学的な意味を説明せよ。  
Explain the scientific meaning of adding the N box in the 4-box model.

### 想定解答 / Expected answers

N box は北大西洋を表し、低緯度から来た水が沈み込む場所を表す。これにより、深層水形成と南大洋湧昇を別々に扱える。  

The N box represents the North Atlantic, where water from low latitudes sinks. This allows deep-water formation and Southern Ocean upwelling to be treated separately.

---

## 課題 3 / Exercise 3

### 問題 / Question

`FHD` を小さくしたとき、大気 pCO2、深層 DIC、深層 O2 はどう変化するか。  
When `FHD` is reduced, what happens to atmospheric pCO2, deep DIC, and deep O2?

### 想定解答 / Expected answers

`FHD` を小さくすると、H と D の交換が弱まり、深層がより孤立する。多くの場合、深層 DIC は増え、深層 O2 は減り、大気 pCO2 は下がる方向に動く。  

Reducing `FHD` weakens exchange between H and D and makes the deep ocean more isolated. In many cases, deep DIC increases, deep O2 decreases, and atmospheric pCO2 tends to decrease.

---

## 課題 4 / Exercise 4

### 問題 / Question

`Tcir` を小さくしたとき、N box の PO4 はどう変化するか。  
When `Tcir` is reduced, what happens to PO4 in the N box?

### 想定解答 / Expected answers

`Tcir` を小さくすると、低緯度から N box へ運ばれる水が減るため、N box の栄養塩供給が変化する。N box は低栄養塩のままになりやすいが、応答は輸出生産やガス交換にも依存する。  

Reducing `Tcir` changes the supply of water from low latitudes to the N box. The N box tends to remain nutrient-poor, but the response also depends on export production and gas exchange.

---

## 課題 5 / Exercise 5

### 問題 / Question

4-box でもまだ表現できない海洋過程を 3 つ挙げよ。  
List three ocean processes that still cannot be represented by the 4-box model.

### 想定解答 / Expected answers

4-box でも、中層水、太平洋と大西洋の違い、インド洋、海氷、季節変化、空間的に連続した南大洋の構造などは表現できない。  

Even the 4-box model cannot represent intermediate waters, Atlantic–Pacific differences, the Indian Ocean, sea ice, seasonality, or the continuous spatial structure of the Southern Ocean.

---

# 04-02 Four-box model

## 課題 1 / Exercise 1

### 問題 / Question

4-box モデルの変数名のルールを説明せよ。  
Explain the naming rule for variables in the 4-box model.

### 想定解答 / Expected answers

変数名は「トレーサー名 + ボックス名」で作られる。例えば `DICN` は N box の DIC、`PO4D` は D box の PO4 を表す。  

Variable names are formed as **tracer name + box name**. For example, `DICN` means DIC in the N box, and `PO4D` means PO4 in the D box.

---

## 課題 2 / Exercise 2

### 問題 / Question

L box の PO4 保存則を言葉で説明せよ。  
Explain the PO4 conservation law for the L box in words.

### 想定解答 / Expected answers

L box の PO4 は、H box からの輸送で増え、L box の生物輸出で減る。  

PO4 in the L box increases through transport from the H box and decreases through biological export in the L box.

---

## 課題 3 / Exercise 3

### 問題 / Question

D box の PO4 保存則では、なぜ `EPH + EPL + EPN` が足されるのか。  
Why is `EPH + EPL + EPN` added in the PO4 conservation equation for the D box?

### 想定解答 / Expected answers

H, L, N の表層で輸出された有機物が深層 D で再無機化されるため、D box では `EPH + EPL + EPN` が足される。  

Organic matter exported from the surface boxes H, L, and N is remineralized in the deep box D, so `EPH + EPL + EPN` is added in the D-box equation.

---

## 課題 4 / Exercise 4

### 問題 / Question

`air_sea=False` にすると、大気 pCO2 と海洋内部のトレーサーはどう変わるか。  
What happens to atmospheric pCO2 and internal ocean tracers when `air_sea=False`?

### 想定解答 / Expected answers

`air_sea=False` では、大気 pCO2 は更新されない。一方、海洋内部では輸送と生物ポンプが残るため、PO4, DIC, O2 の分布は変化しうる。  

When `air_sea=False`, atmospheric pCO2 is not updated. However, internal ocean tracers such as PO4, DIC, and O2 can still change due to transport and the biological pump.

---

## 課題 5 / Exercise 5

### 問題 / Question

生物ポンプを止めると、PO4 と DIC の鉛直差はどう変わるか。  
What happens to the vertical gradients of PO4 and DIC when the biological pump is turned off?

### 想定解答 / Expected answers

生物ポンプを止めると、表層から深層への物質輸送が弱くなるため、PO4 と DIC の鉛直差は小さくなる。  

When the biological pump is turned off, transport of material from the surface to the deep ocean weakens, so the vertical gradients of PO4 and DIC become smaller.

---

# 04-03 Four-box model

## 課題 1 / Exercise 1

### 問題 / Question

`Tcir` を小さくしたとき、N box の PO4 と深層 O2 はどう変化したか。  
When `Tcir` was reduced, how did PO4 in the N box and deep O2 change?

### 想定解答 / Expected answers

`Tcir` を小さくすると、N box と深層 D のつながりが弱くなり、N box の栄養塩供給や深層 O2 が変化する。応答は生物ポンプとガス交換にも依存する。  

Reducing `Tcir` weakens the connection between the N box and the deep box D, affecting nutrient supply to the N box and deep O2. The response also depends on biological export and air–sea gas exchange.

---

## 課題 2 / Exercise 2

### 問題 / Question

`FHD` を小さくしたとき、大気 pCO2、深層 DIC、深層 O2 はどう変化したか。  
When `FHD` was reduced, how did atmospheric pCO2, deep DIC, and deep O2 change?

### 想定解答 / Expected answers

`FHD` を小さくすると、H と D の交換が弱まり、深層が孤立しやすくなる。典型的には深層 DIC が増え、深層 O2 が下がり、大気 pCO2 は下がる方向に動く。  

Reducing `FHD` weakens exchange between H and D and makes the deep ocean more isolated. Typically, deep DIC increases, deep O2 decreases, and atmospheric pCO2 tends to decrease.

---

## 課題 3 / Exercise 3

### 問題 / Question

`LF` と `CEPH` の効果は同じか。違う場合、なぜ違うのか。  
Are the effects of `LF` and `CEPH` the same? If not, why are they different?

### 想定解答 / Expected answers

`LF` は低緯度 L の生物ポンプ、`CEPH` は南大洋 H の生物ポンプを変える。H は深層と大気の接点なので、`CEPH` は深層水由来の CO2 放出に効きやすい。一方、`LF` は広い低緯度表層での DIC 除去に効く。  

`LF` changes the biological pump in the low-latitude box L, whereas `CEPH` changes the biological pump in the Southern Ocean box H. Because H is the contact point between deep water and the atmosphere, `CEPH` strongly affects CO2 outgassing from deep water. In contrast, `LF` mainly affects DIC removal in the broad low-latitude surface ocean.

---

## 課題 4 / Exercise 4

### 問題 / Question

`RRC` を変えると、なぜ DIC だけでなく ALK や CO3 も変わるのか。  
Why do ALK and CO3, not only DIC, change when `RRC` is varied?

### 想定解答 / Expected answers

`RRC` は CaCO3 輸出の比を表すため、DIC だけでなく ALK も変える。ALK が変わると炭酸イオン濃度や pCO2 も変化する。  

`RRC` represents the ratio of CaCO3 export, so it changes not only DIC but also ALK. Changes in ALK also affect carbonate ion concentration and pCO2.

---

## 課題 5 / Exercise 5

### 問題 / Question

大気 pCO2 を下げるために最も有効そうなパラメタを 1 つ選び、理由を説明せよ。  
Choose one parameter that seems most effective for lowering atmospheric pCO2 and explain why.

### 想定解答 / Expected answers

このモデルでは、`FHD` や `CEPH` が大気 pCO2 に効きやすい。`FHD` を弱めると深層炭素が大気へ出にくくなり、`CEPH` を強めると南大洋表層から炭素を深層へ輸送しやすくなるためである。  

In this model, `FHD` and `CEPH` are among the most effective parameters controlling atmospheric pCO2. Weakening `FHD` makes it more difficult for deep-ocean carbon to reach the atmosphere, whereas increasing `CEPH` enhances carbon transport from the Southern Ocean surface to the deep ocean.

---

# 04-04 Four-box model

## 課題 1 / Exercise 1

### 問題 / Question

なぜモデルを改造するとき、元の関数を直接書き換えず、コピーを作るべきなのか。  
Why should we copy a function instead of directly modifying the original function?

### 想定解答 / Expected answers

元の関数を直接書き換えると、標準実験を再現できなくなる。コピーを作れば、標準実験と改造実験を比較でき、バグが出たときにも戻りやすい。  
If the original function is directly modified, the standard experiment may no longer be reproducible. By making a copy, we can compare the standard and modified experiments and recover more easily from bugs.

---

## 課題 2 / Exercise 2

### 問題 / Question

時間変化する `FHD` を入れることで、どのような科学的問題を表現できるか。  
What scientific problem can be represented by introducing time-varying `FHD`?

### 想定解答 / Expected answers

退氷期や氷期のように、南大洋の深層ベンチレーションが時間とともに変化する問題を表現できる。  
It can represent problems such as glacial or deglacial changes in which Southern Ocean ventilation changes through time.

---

## 課題 3 / Exercise 3

### 問題 / Question

生物ポンプに温度依存性を入れると、どの box の輸出生産が特に変わりやすいか。  
When temperature dependence is added to the biological pump, which box's export production changes most easily?

### 想定解答 / Expected answers

温度依存性を入れると、暖かい低緯度 L box の輸出生産が特に強く変化しやすい。  
With temperature dependence, export production in the warm low-latitude L box tends to change strongly.

---

## 課題 4 / Exercise 4

### 問題 / Question

南大洋 H のガス交換を弱める実験は、どのような現実過程の簡略表現と考えられるか。  
Weakening gas exchange in H can be interpreted as a simplified representation of what real-world process?

### 想定解答 / Expected answers

南大洋 H のガス交換を弱めることは、海氷の拡大、風速低下、成層化による大気海洋交換の抑制などの簡略表現と考えられる。  
Weakening gas exchange in H can be viewed as a simplified representation of sea-ice expansion, weaker winds, or stronger stratification suppressing air-sea exchange.

---

## 課題 5 / Exercise 5

### 問題 / Question

自分のシナリオ実験について、問い・仮説・結果・解釈を 5〜10 行でまとめよ。  
Summarize your own scenario experiment in 5–10 lines: question, hypothesis, result, and interpretation.

### 想定解答 / Expected answers

南大洋ベンチレーションを弱め、生物ポンプを強めると、大気 CO2 は低下した。これは、深層に蓄積した炭素が大気へ出にくくなり、さらに表層炭素が深層へ輸送されたためと解釈できる。ただし、このモデルは空間構造や海氷、鉄制限を明示的に含まないため、現実の氷期 CO2 変化を完全に説明するものではない。  
When Southern Ocean ventilation was weakened and the biological pump was strengthened, atmospheric CO2 decreased. This can be interpreted as deep carbon becoming less able to reach the atmosphere, while surface carbon was additionally exported to the deep ocean. However, this model does not explicitly include spatial structure, sea ice, or iron limitation, so it cannot fully explain real glacial CO2 change.

---

# 04-05 Four-box model

## 課題 1 / Exercise 1

### 問題 / Question

3-box から 4-box へ拡張することで、何が新しく表現できるようになったか。  
What becomes newly representable by extending from 3-box to 4-box?

### 想定解答 / Expected answers

4-box では、北大西洋の沈み込みと南大洋の湧昇を別々に表現できる。これにより、深層水形成と南大洋ベンチレーションが大気 CO2 に与える影響を分けて調べられる。  
The 4-box model can separately represent North Atlantic sinking and Southern Ocean upwelling. This allows us to examine the effects of deep-water formation and Southern Ocean ventilation on atmospheric CO2 separately.

---

## 課題 2 / Exercise 2

### 問題 / Question

H box と N box の役割の違いを説明せよ。  
Explain the difference between the roles of the H box and the N box.

### 想定解答 / Expected answers

H box は南大洋表層を表し、深層水の湧昇と大気海洋 CO2 交換に関係する。N box は北大西洋表層を表し、低緯度から来た水が沈み込む場所を表す。  
The H box represents the Southern Ocean surface and is related to deep-water upwelling and air-sea CO2 exchange. The N box represents the North Atlantic surface, where water from low latitudes sinks.

---

## 課題 3 / Exercise 3

### 問題 / Question

`FHD` を弱めると、大気 pCO2 が低下しやすい理由を説明せよ。  
Explain why weakening `FHD` tends to lower atmospheric pCO2.

### 想定解答 / Expected answers

`FHD` を弱めると、深層 D と南大洋 H の交換が弱くなる。そのため、深層に蓄積した DIC が南大洋表層へ上がりにくくなり、大気へ CO2 が放出されにくくなる。  
Weakening `FHD` reduces exchange between deep D and Southern Ocean H. Therefore, DIC stored in the deep ocean is less able to reach the Southern Ocean surface and outgas to the atmosphere.

---

## 課題 4 / Exercise 4

### 問題 / Question

生物ポンプを強めると、大気 pCO2 が低下しやすい理由を説明せよ。  
Explain why strengthening the biological pump tends to lower atmospheric pCO2.

### 想定解答 / Expected answers

生物ポンプを強めると、表層の DIC と栄養塩が有機物として取り除かれ、深層へ輸送される。表層 pCO2 が低くなりやすいため、大気 pCO2 も低下しやすい。  
A stronger biological pump removes DIC and nutrients from the surface as organic matter and transfers them to the deep ocean. Surface pCO2 tends to become lower, so atmospheric pCO2 also tends to decrease.

---

## 課題 5 / Exercise 5

### 問題 / Question

4-box モデルの限界を 3 つ挙げよ。  
List three limitations of the 4-box model.

### 想定解答 / Expected answers

4-box モデルでは、太平洋とインド洋の違い、中層水、季節変化、海氷、風の空間構造、鉄制限、生態系の多様性などは明示的に表現できない。  
The 4-box model cannot explicitly represent Pacific-Indian differences, intermediate waters, seasonality, sea ice, spatial wind structure, iron limitation, or ecosystem diversity.

---

## 課題 6 / Exercise 6

### 問題 / Question

短い研究報告として、標準実験と複合シナリオの違いを 5〜10 行で説明せよ。  
As a short research report, explain the difference between the standard experiment and the combined scenario in 5–10 lines.

### 想定解答 / Expected answers

解答例: 複合シナリオでは、南大洋ベンチレーションの弱化、生物ポンプの強化、南大洋ガス交換の弱化を同時に与えた。その結果、標準実験よりも大気 pCO2 が大きく低下した。これは、深層炭素が大気へ戻りにくくなり、さらに表層炭素が深層へ輸送されやすくなったためと解釈できる。ただし、このモデルは単純化されており、現実の氷期 CO2 低下を完全に再現するものではない。  
Example: In the combined scenario, Southern Ocean ventilation was weakened, the biological pump was strengthened, and Southern Ocean gas exchange was reduced simultaneously. As a result, atmospheric pCO2 decreased more strongly than in the standard experiment. This can be interpreted as deep carbon becoming less able to return to the atmosphere, while surface carbon was more efficiently transferred to the deep ocean. However, this model is simplified and does not fully reproduce real glacial CO2 decrease.

---

# xx-x TBD

## 課題 1 / Exercise 1

### 問題 / Question

### 想定解答 / Expected answers

---

## 課題 2 / Exercise 2

### 問題 / Question

### 想定解答 / Expected answers

---

## 課題 3 / Exercise 3

### 問題 / Question

### 想定解答 / Expected answers

---

## 課題 4 / Exercise 4

### 問題 / Question

### 想定解答 / Expected answers

---

## 課題 5 / Exercise 5

### 問題 / Question

### 想定解答 / Expected answers

---

### 課題 1 / Exercise 1

### 課題 2 / Exercise 2

### 課題 3 / Exercise 3

### 課題 4 / Exercise 4

### 課題 5 / Exercise 5

### 課題 6 / Exercise 6

---

## 想定解答 / Expected answers

### 課題 1

### 課題 2

### 課題 3

### 課題 4

### 課題 5

### 課題 6
]
