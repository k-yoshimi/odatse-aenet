チュートリアル: AENET 学習パイプライン
======================================

本チュートリアルでは、N2 二量体を例に、第一原理計算（Quantum ESPRESSO）から
教師データを作成し、AENET で機械学習ポテンシャルを構築する一連の手順を説明します。

サンプルファイルは ``sample/aenet_training/`` にあります。

ワークフロー概要
----------------

.. code-block:: text

   1. 構造生成 (relax/)
      └→ ランダムなN2二量体構造を生成
      └→ QE で構造緩和・エネルギー計算

   2. 教師データ作成 (relax/)
      └→ QE 出力からエネルギー・力を抽出
      └→ XSF 形式の教師データを生成

   3. フィンガープリント生成 (generate/)
      └→ AENET generate.x で原子記述子を計算

   4. ニューラルネットワーク学習 (train/)
      └→ AENET train.x でポテンシャルを学習

   5. 予測・検証 (predict/)
      └→ AENET predict.x でテスト構造のエネルギーを予測
      └→ QE 計算との比較

Step 1: 構造生成
-----------------

``sample/aenet_training/relax/`` ディレクトリで作業します。

``structure_make.py`` は、N-N 結合距離を 0.5〜2.0 Å の範囲でランダムに生成し、
Quantum ESPRESSO の入力ファイルを作成します。

.. code-block:: python

   import random
   structure_list = []
   for i in range(20):
       structure_list.append(random.uniform(0.5, 2.0))

テンプレートファイル (``template.txt``) は QE の入力形式で、``value_01`` が
N-N 結合距離のプレースホルダーです:

.. code-block:: text

   &CONTROL
     calculation = 'relax'
     pseudo_dir  = './'
     outdir      = './'
   /
   &SYSTEM
     ntyp = 1, nat = 2, ibrav = 0
     ecutwfc = 44, ecutrho = 320
     occupations = 'smearing'
     smearing    = 'mp'
     degauss     = 0.01
   /
   ...
   ATOMIC_POSITIONS angstrom
     N 0.00 0.00 0.00
     N 0.00 0.00 value_01

実行:

.. code-block:: bash

   python3 structure_make.py

これにより ``directory_0/`` 〜 ``directory_19/`` が生成され、それぞれに QE 入力ファイルが配置されます。

各ディレクトリで QE を実行します:

.. code-block:: bash

   srun pw.x < n2_dimer.pwi > n2_dimer.pwo

Step 2: 教師データ作成
-----------------------

QE の出力ファイルからエネルギーと原子間力を抽出し、AENET の XSF 形式に変換します。

``teach_data_make.py`` は以下を行います:

- ``.pwo`` ファイルからエネルギーを抽出（Ry → eV 変換: ``ase.units.Ry``）
- 原子座標と力を抽出（力の変換: ``ase.units.Ry / 0.529177`` [Ry/Bohr → eV/Å]）
- XSF 形式で出力

.. code-block:: bash

   cd directory_0
   python3 ../teach_data_make.py --input n2_dimer.pwo --output-dir teach_data

出力例（teach_data_1.xsf）:

.. code-block:: text

   # total energy = -767.51 eV

   ATOMS
   N   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000  -12.34567890
   N   0.00000000   0.00000000   1.04970000   0.00000000   0.00000000   12.34567890

Step 3: フィンガープリント生成
-------------------------------

``sample/aenet_training/generate/`` ディレクトリで作業します。

``generate.in`` で教師データのパスと原子種を指定します:

.. code-block:: text

   OUTPUT N2.train

   TYPES
   1
   N 389.83

   SETUPS
   N N.fingerprint.stp

   FILES
   20
   ../relax/directory_0/teach_data/teach_data_1.xsf
   ../relax/directory_1/teach_data/teach_data_1.xsf
   ...

実行:

.. code-block:: bash

   generate.x generate.in > generate.out

出力:

- ``N2.train``: 学習用バイナリデータセット
- ``N.fingerprint.stp``: フィンガープリント設定ファイル

Step 4: ニューラルネットワーク学習
-----------------------------------

``sample/aenet_training/train/`` ディレクトリで作業します。

``train.in`` で学習パラメータを設定します:

.. code-block:: text

   TRAININGSET N2.train
   TESTPERCENT 10
   ITERATIONS  3000

   bfgs

   NETWORKS
   # atom   network           hidden
   # types  file-name         layers   nodes:activation
   N        N.5t-5t.ann       2        10:tanh 10:tanh

主な設定:

- **TESTPERCENT 10**: 教師データの 10% をテストデータとして使用
- **ITERATIONS 3000**: BFGS 最適化を 3000 回反復
- **ネットワーク構造**: 入力層 → 10 ノード(tanh) → 10 ノード(tanh) → 出力層

実行:

.. code-block:: bash

   train.x train.in > train.out

学習結果の例（train.out より）:

.. code-block:: text

   Number of training structures :         18
   Number of testing structures  :          2
   Total number of weights       :        221
   Atomic energy shift           :   -747.530391 eV

出力:

- ``N.5t-5t.ann``: 学習済み ANN ポテンシャルファイル

Step 5: 予測・検証
-------------------

``sample/aenet_training/predict/`` ディレクトリで作業します。

テスト構造（N-N 距離 0.00〜2.00 Å、0.01 Å 刻みの 201 構造）に対して
学習済みポテンシャルによるエネルギー予測を行います。

.. code-block:: bash

   predict.x predict.in > predict.out

結果の可視化
~~~~~~~~~~~~

``plot_distance_energy.py`` を使用して、ANN ポテンシャルによる予測エネルギーと QE の
第一原理計算結果を比較するプロットを作成できます:

.. code-block:: bash

   python3 plot_distance_energy.py \
       --predict-out predict.out \
       --train-out ../train/train.out \
       --qe-dir ../relax \
       --n-structures 20 \
       --output distance_energy_plot.png

このスクリプトは以下の 3 つのデータセットをプロットします:

- **ANN 予測** (青): ``predict.out`` からの原子間距離 vs エネルギー
- **QE 教師データ** (赤): 緩和計算の XSF ファイルからの結果
- **QE テストデータ** (緑): ``train.out`` のテストセットに含まれるデータ点

計算結果
~~~~~~~~

以下のグラフは、ANN ポテンシャルによる予測エネルギーと QE の第一原理計算結果を比較したものです。

.. figure:: ../../common/img/QE_ANN.png
   :width: 80%
   :align: center

   N2 二量体のエネルギー-距離曲線。青: ANN 予測、赤: QE 教師データ、緑: QE テストデータ。

主な結果:

- **ANN 最安定距離**: 約 1.11 Å（エネルギー: -767.81 eV）
- **QE 最安定距離**: 約 1.05 Å（エネルギー: -767.51 eV）
- 教師データ（20 点）の範囲内で ANN ポテンシャルが QE の結果をよく再現
