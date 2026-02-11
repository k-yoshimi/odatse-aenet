チュートリアル: EMT による Ni バルク分子動力学
===============================================

本チュートリアルでは、ASE に内蔵された EMT (Effective Medium Theory) ポテンシャルを用いて
Ni バルクの NVT 分子動力学シミュレーションを実行する方法を説明します。

外部の計算エンジン（QE や ELSES）を必要とせず、ASE だけで完結するため、
MD シミュレーションの基本的なワークフローを確認するのに適しています。

サンプルファイルは ``sample/emt_md/`` にあります。

前提条件
--------

- ASE がインストール済みであること

ファイル構成
------------

.. code-block:: text

   sample/emt_md/
   └── run_md.py   # MD シミュレーションスクリプト

シミュレーション設定
--------------------

``run_md.py`` は 2 段階の温度制御で NVT MD を実行します:

1. **平衡化**: 低温（10 K）で 200 ステップ
2. **本計算**: 高温（500 K）に昇温して 400 ステップ

構造の準備
~~~~~~~~~~

.. code-block:: python

   from ase.build import bulk, make_supercell

   atoms = bulk("Ni", "fcc", a=3.5, cubic=True)
   atoms = make_supercell(atoms, np.diag([3, 3, 3]))
   atoms.calc = EMT()

FCC Ni の単位胞を 3×3×3 のスーパーセル（108 原子）に拡張し、EMT ポテンシャルを設定します。

MD 設定
~~~~~~~

.. code-block:: python

   from ase import units
   from ase.md.nvtberendsen import NVTBerendsen
   from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

   dt = 2.0 * units.fs
   MaxwellBoltzmannDistribution(atoms, temperature_K=10)

   dyn = NVTBerendsen(atoms, dt, temperature=10, taut=20*units.fs,
                      trajectory="md.traj")
   dyn.run(200)       # 平衡化

   dyn.set_temperature(500)
   dyn.run(400)       # 本計算

- **アンサンブル**: NVT（Berendsen サーモスタット）
- **タイムステップ**: 2 fs
- **温度制御時定数**: 20 fs
- **初期温度**: 10 K → 500 K に昇温

実行方法
--------

.. code-block:: bash

   cd sample/emt_md
   python3 run_md.py

コマンドラインオプションで各パラメータを変更できます:

.. code-block:: bash

   python3 run_md.py --element Ni --supercell 3 \
       --temp-equil 10 --nsteps-equil 200 \
       --temp-prod 500 --nsteps-prod 400

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - オプション
     - デフォルト値
     - 説明
   * - ``--element``
     - ``Ni``
     - 元素記号
   * - ``--supercell``
     - ``3``
     - スーパーセルの倍率（各方向同一）
   * - ``--temp-equil``
     - ``10``
     - 平衡化温度 [K]
   * - ``--nsteps-equil``
     - ``200``
     - 平衡化ステップ数
   * - ``--temp-prod``
     - ``500``
     - 本計算温度 [K]
   * - ``--nsteps-prod``
     - ``400``
     - 本計算ステップ数
   * - ``--traj``
     - ``md.traj``
     - 出力トラジェクトリファイル

出力
----

- 各ステップで時間と温度がコンソールに出力されます
- ``md.traj``: ASE トラジェクトリファイル

.. code-block:: text

   --- Equilibration at 10 K ---
   time=    0 fs T= 10 K
   time=   40 fs T= 12 K
   ...
   --- Production at 500 K ---
   time=  400 fs T= 489 K
   time=  440 fs T= 503 K
   ...

平衡化段階で 10 K 付近に安定化した後、500 K に昇温すると Berendsen サーモスタットにより
目標温度に緩和していく様子が確認できます。

.. note::

   EMT は Cu, Ag, Au, Ni, Pd, Pt, C, N, O に対応しています。
   他の元素を使用する場合は、適切な計算エンジン（QE 等）に切り替えてください。
