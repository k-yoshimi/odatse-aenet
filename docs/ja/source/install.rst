インストール
============

必要環境
--------

- Python >= 3.9
- numpy >= 1.22
- `ODAT-SE <https://github.com/issp-center-dev/ODAT-SE>`_ >= 3.0
- `AENET <http://ann.atomistic.net/>`_ (predict.x, generate.x, train.x)

サンプルによっては以下も必要です:

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ (Atomic Simulation Environment)
- `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_ (pw.x)
- `ELSES <https://www.elses.jp/>`_

ODAT-SE のインストール
-----------------------

PyPI からインストールする場合:

.. code-block:: bash

   pip3 install odat-se[all]

ソースからインストールする場合:

.. code-block:: bash

   git clone https://github.com/issp-center-dev/ODAT-SE.git
   cd ODAT-SE
   pip3 install .[all]

AENET のインストール
---------------------

AENET（Atomic Energy Network）は機械学習ポテンシャルの構築に使用します。
ビルドには Fortran コンパイラ（gfortran）と LAPACK/BLAS が必要です。

macOS の場合の前準備
~~~~~~~~~~~~~~~~~~~~

gfortran が入っていない場合は Homebrew でインストールします:

.. code-block:: bash

   brew install gcc

これにより ``gfortran`` コマンドが利用可能になります。
macOS では LAPACK/BLAS として Accelerate フレームワークが標準で利用できます。

ソースコードの取得
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/atomisticnet/aenet.git
   cd aenet

L-BFGS-B ライブラリのビルド
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AENET 本体のビルド前に、同梱されている L-BFGS-B ライブラリを先にビルドする必要があります:

.. code-block:: bash

   cd lib
   make static
   cd ..

``lib/liblbfgsb.a`` が生成されていることを確認してください。

AENET 本体のビルド
~~~~~~~~~~~~~~~~~~

``src/`` ディレクトリで、環境に合った Makefile を指定してビルドします:

.. code-block:: bash

   cd src
   make -f makefiles/Makefile.gfortran_serial_MacOS   # macOS (Apple Silicon / Intel)
   cd ..

.. note::

   Linux 環境の場合は、MPI の有無に応じて以下を選択してください:

   - ``Makefile.gfortran_serial`` — シリアル版
   - ``Makefile.gfortran_mpi`` — MPI 並列版
   - ``Makefile.gfortran_openblas_serial`` — OpenBLAS 利用版

   利用可能な Makefile 一覧は ``src/makefiles/`` を参照してください。

ビルドが完了すると ``bin/`` 以下に以下の実行ファイルが生成されます:

- ``generate.x`` — 学習データ（構造記述子）の生成
- ``train.x`` — ニューラルネットワークポテンシャルの学習
- ``predict.x`` — 学習済みポテンシャルによるエネルギー予測

.. note::

   実行ファイル名にはバージョンとコンパイラのサフィックスが付きます
   （例: ``predict.x-2.0.4-gfortran_serial``）。
   必要に応じてシンボリックリンクを作成してください:

   .. code-block:: bash

      cd bin
      ln -s predict.x-2.0.4-gfortran_serial predict.x
      ln -s generate.x-2.0.4-gfortran_serial generate.x
      ln -s train.x-2.0.4-gfortran_serial train.x

パスの設定:

.. code-block:: bash

   export PATH=/path/to/aenet/bin:$PATH

動作確認:

.. code-block:: bash

   generate.x
   # "generate.x - training set generation" と表示されれば成功

.. note::

   詳細は `AENET 公式サイト <http://ann.atomistic.net/>`_ を参照してください。

ASE のインストール
-------------------

ASE（Atomic Simulation Environment）は分子動力学シミュレーションや
各種計算エンジンとの連携に使用します。

.. code-block:: bash

   pip3 install ase

インストール確認:

.. code-block:: bash

   python3 -c "import ase; print(ase.__version__)"

.. note::

   詳細は `ASE 公式ドキュメント <https://wiki.fysik.dtu.dk/ase/>`_ を参照してください。

Quantum ESPRESSO のインストール
---------------------------------

Quantum ESPRESSO は第一原理計算エンジンで、フラーレンMDチュートリアルで使用します。

前準備（macOS）
~~~~~~~~~~~~~~~

cmake が必要です:

.. code-block:: bash

   brew install cmake

ソースからビルド
~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone --depth 1 --branch qe-7.4 https://github.com/QEF/q-e.git qe-7.4
   cd qe-7.4

.. warning::

   v7.3.1 以前は GCC 15 (gfortran 15) で mbd ライブラリのビルドに失敗します。
   **v7.4 以降の使用を推奨します。**

cmake でビルドを構成します。MPI なし・FFTW 内蔵版の場合:

.. code-block:: bash

   cmake -B build \
     -DCMAKE_Fortran_COMPILER=gfortran \
     -DCMAKE_C_COMPILER=gcc-15 \
     -DQE_ENABLE_MPI=OFF \
     -DQE_FFTW_VENDOR=Internal

.. note::

   - ``gcc-15`` の部分は環境に合わせて変更してください
     （例: Linux で ``gcc`` が GNU の場合はそのまま ``gcc``）。
   - MPI 並列版が必要な場合は ``-DQE_ENABLE_MPI=OFF`` を外し、
     ``-DCMAKE_Fortran_COMPILER=mpif90`` を指定してください。
   - 外部 FFTW3 がインストール済みの場合は ``-DQE_FFTW_VENDOR=Internal`` は不要です。

pw.x のみをビルドする場合:

.. code-block:: bash

   cmake --build build --target pw -j4

全体をビルドする場合:

.. code-block:: bash

   cmake --build build -j4

ビルド後、``pw.x`` にパスを通します:

.. code-block:: bash

   export PATH=/path/to/qe-7.4/build/bin:$PATH

動作確認:

.. code-block:: bash

   pw.x --version
   # "Program PWSCF v.7.4" と表示されれば成功

.. note::

   - 擬ポテンシャルファイル（``.UPF``）が必要です。
     `SSSP ライブラリ <https://www.materialscloud.org/discover/sssp/>`_
     などから取得できます。
   - 詳細は `Quantum ESPRESSO 公式サイト <https://www.quantum-espresso.org/>`_ を参照してください。

ELSES のインストール
---------------------

ELSES は大規模電子状態計算プログラムで、フラーレンMDチュートリアル（ELSES版）で使用します。

ELSES のソースコードは `ELSES コンソーシアム <http://www.elses.jp/>`_ のメンバーに提供されています。
利用するには、コンソーシアムへの参加が必要です。

詳細は `ELSES 公式サイト <http://www.elses.jp/>`_ を参照してください。

odatse-aenet のインストール
----------------------------

.. code-block:: bash

   git clone https://github.com/k-yoshimi/odatse-aenet.git
   cd odatse-aenet
   pip3 install .

インストール後、``odatse-aenet`` コマンドが利用可能になります。

動作確認
--------

.. code-block:: bash

   odatse-aenet --version

バージョン番号が表示されれば、正しくインストールされています。
