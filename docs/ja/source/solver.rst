AENET ソルバー
==============

概要
----

AenetSolver は、ODAT-SE の ``odatse.solver.SolverBase`` を継承したソルバークラスです。
AENET の ``predict.x`` をサブプロセスとして実行し、構造テンプレートにパラメータを代入して
原子間エネルギーを計算します。

ソルバーの仕組み
----------------

``evaluate`` メソッドが呼ばれると、以下の処理が行われます:

1. 作業ディレクトリ ``Log{step}_{iset}/`` を作成
2. ``predict.in`` と ANN ポテンシャルファイルをコピー
3. テンプレートファイル内のプレースホルダー（例: ``value_01``）をパラメータ値で置換し、入力構造ファイルを生成
4. ``predict.x`` を実行
5. ``predict.out`` から全エネルギーと原子数を抽出
6. エネルギー/原子数（eV/atom）を R-factor として返却

ソースコード構造
~~~~~~~~~~~~~~~~

.. code-block:: text

   src/AenetSolver/
     __init__.py        # Solver クラスのエクスポート
     _main.py           # CLI エントリポイント (odatse.initialize + algorithm dispatch)
     aenet.py           # Solver クラス本体

入力ファイル
------------

input.toml
~~~~~~~~~~

ODAT-SE 標準の TOML 設定ファイルです。

.. code-block:: toml

   [base]
   dimension = 1
   output_dir = "output"

   [solver]
   name = "aenet"

   [solver.config]
   aenet_exec_file = "predict.x"           # predict.x のパス
   aenet_ann_potential = "N.5t-5t.ann"      # ANN ポテンシャルファイル
   aenet_template_file = "template.xsf"    # テンプレート構造ファイル (デフォルト)
   aenet_input_file = "structure_data.xsf"  # 生成される入力ファイル名 (デフォルト)
   bulk_output_file = "predict.in"          # AENET 設定ファイル (デフォルト)

   [solver.param]
   string_list = ["value_01"]              # テンプレート内のプレースホルダー

   [solver.post]
   remove_work_dir = false                 # 計算後に作業ディレクトリを削除するか

   [algorithm]
   name = "mapper"                          # 使用するアルゴリズム
   label_list = ["z"]

   [algorithm.param]
   min_list = [0.8]
   max_list = [1.4]
   num_list = [101]

solver.config セクション
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - パラメータ
     - デフォルト値
     - 説明
   * - ``aenet_exec_file``
     - ``predict.x``
     - AENET predict.x 実行ファイルのパス
   * - ``aenet_ann_potential``
     - ``N.10t-10t.ann``
     - 学習済み ANN ポテンシャルファイル
   * - ``aenet_template_file``
     - ``template.xsf``
     - 構造テンプレートファイル
   * - ``aenet_input_file``
     - ``structure_data.xsf``
     - 生成される構造入力ファイル名
   * - ``bulk_output_file``
     - ``predict.in``
     - AENET 設定ファイル

solver.param セクション
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - パラメータ
     - デフォルト値
     - 説明
   * - ``string_list``
     - ``["value_01", "value_02"]``
     - テンプレート内で置換するプレースホルダーのリスト

テンプレートファイル
~~~~~~~~~~~~~~~~~~~~

テンプレートファイル (``template.xsf``) は XSF 形式の構造ファイルで、
最適化パラメータの位置にプレースホルダーを記述します。

.. code-block:: text

   ATOMS
   N             0.0000000000        0.0000000000        0.0000000000
   N             0.0000000000        0.0000000000        value_01

この例では ``value_01`` が N-N 結合距離に対応し、アルゴリズムによって最適化されます。

predict.in
~~~~~~~~~~

AENET 予測計算の設定ファイルです。

.. code-block:: text

   TYPES
   1
   N

   NETWORKS
     N N.5t-5t.ann

   FORCES

   FILES
   1
   structure_data.xsf

実行方法
--------

CLI 実行:

.. code-block:: bash

   odatse-aenet input.toml

MPI 並列実行:

.. code-block:: bash

   mpiexec -np 4 odatse-aenet input.toml

スクリプトから直接実行:

.. code-block:: bash

   python3 src/main.py input.toml

