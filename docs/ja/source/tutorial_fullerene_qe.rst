チュートリアル: フラーレン分子動力学（Quantum ESPRESSO）
=========================================================

本チュートリアルでは、ASE (Atomic Simulation Environment) と
Quantum ESPRESSO を用いて、フラーレン（C60 系）の NVT 分子動力学シミュレーションを
実行する方法を説明します。

サンプルファイルは ``sample/fullerene_qe_md/`` にあります。

前提条件
--------

- Quantum ESPRESSO (pw.x) がインストール済みであること
- ASE がインストール済みであること
- C の擬ポテンシャルファイル（``C.pbe-n-kjpaw_psl.1.0.0.UPF``）
- 結晶構造ファイル（``C.cif``）

ファイル構成
------------

.. code-block:: text

   sample/fullerene_qe_md/
   ├── do.sh                 # 実行スクリプト（ローカル実行用）
   ├── job_ohtaka.sh         # ISSP Ohtaka クラスタ用ジョブスクリプト
   ├── ase_qe.py             # MD シミュレーションスクリプト
   └── time_energy_plot.py   # 結果可視化スクリプト

シミュレーション設定
--------------------

``ase_qe.py`` の主要な設定項目:

構造の読み込みとスーパーセル作成
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ase.io import read
   from ase.build import make_supercell

   atoms = read("C.cif")
   atoms.pbc = [1, 1, 1]
   supercell_matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
   supercell = make_supercell(atoms, supercell_matrix)

CIF ファイルから構造を読み込み、x 方向に 2 倍のスーパーセル（480 原子）を作成します。

Quantum ESPRESSO の設定
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   input_data = {
       'control': {
           'calculation': 'scf',
           'pseudo_dir': './',
           'outdir': './'
       },
       'system': {
           'ecutwfc': 40,    # 波動関数カットオフ [Ry]
           'ecutrho': 160,   # 電荷密度カットオフ [Ry]
           'occupations': 'smearing',
           'smearing': 'gaussian',
           'degauss': 0.02,
       },
       'electrons': {
           'conv_thr': 1e-7
       }
   }

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - パラメータ
     - 説明
   * - ``ecutwfc``
     - 波動関数のカットオフエネルギー (40 Ry)
   * - ``ecutrho``
     - 電荷密度のカットオフエネルギー (160 Ry)
   * - ``smearing``
     - ガウシアンスメアリング (degauss = 0.02 Ry)
   * - ``conv_thr``
     - SCF 収束閾値 (10\ :sup:`-7` Ry)

MD 設定
~~~~~~~

.. code-block:: python

   from ase import units
   from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
   from ase.md.nvtberendsen import NVTBerendsen

   dt = 2 * units.fs          # タイムステップ: 2 fs
   temperature_K = 300        # 目標温度: 300 K
   nsteps = 100               # ステップ数: 100

   MaxwellBoltzmannDistribution(supercell, temperature_K=300)

   taut = 0.5 * 10 * units.fs  # 温度制御時定数: 5 fs
   dyn = NVTBerendsen(supercell, timestep=dt, temperature=temperature_K,
                      taut=taut, logfile='md.log')
   dyn.run(nsteps)

- **アンサンブル**: NVT（Berendsen サーモスタット）
- **初期速度**: マクスウェル・ボルツマン分布（300 K）
- **タイムステップ**: 2 fs
- **総シミュレーション時間**: 200 fs (100 ステップ × 2 fs)

実行方法
--------

ローカル実行:

.. code-block:: bash

   cd sample/fullerene_qe_md
   sh do.sh

ISSP Ohtaka クラスタ上で実行する場合:

.. code-block:: bash

   sbatch job_ohtaka.sh

``job_ohtaka.sh`` はパーティション、ノード数、QE モジュールの読み込み等を含む
SLURM ジョブスクリプトです。環境に合わせて編集してください。

出力ファイル:

- ``md.log``: 各ステップの時間、エネルギー、運動エネルギー、温度
- ``benzene_optimization.traj``: ASE トラジェクトリファイル
- ``espresso.pwi`` / ``espresso.pwo``: QE の入出力

結果の可視化
------------

``time_energy_plot.py`` で運動エネルギーの時間変化をプロットします:

.. code-block:: bash

   python3 time_energy_plot.py

.. figure:: ../../common/img/time_energy_plot_QE.png
   :width: 80%
   :align: center

   フラーレンスーパーセルの運動エネルギー（eV/atom）の時間変化。
   赤破線は 300 K における理論値 :math:`E_{kin} = \frac{3}{2} k_B T` = 0.0388 eV/atom。

計算結果
~~~~~~~~

md.log の出力形式:

.. code-block:: text

   Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]
   0.0000         -XXX.XX      -XXX.XX       XX.XX     300.0
   0.0020         -XXX.XX      -XXX.XX       XX.XX     2XX.X
   ...

- シミュレーション初期は温度が揺らぎながら 300 K 付近に収束
- Berendsen サーモスタットにより目標温度に緩和
- 十分なステップ数で熱平衡に到達

教師データの生成
----------------

MD トラジェクトリから AENET の教師データを生成するには、
:doc:`tutorial_fullerene_elses` の ``teach_data_make.py`` と同様の手順で
XSF 形式に変換できます。
