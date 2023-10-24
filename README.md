# README

大家好，我是這個 project 的作者，我叫 **謝岳璋**，指導教授是 **劉一宇** 博士。
這個 project 是我的論文 "[Topological Escape Routing for Package Substrate Design Planning](./documentation/_2022__Yueh_chang_Shieh__Topological_Escape_Routing_for_Package_Substrate_Design_Planning.pdf)" 的實作，詳細內容可以參考我的論文。

![The Industrial 2 Complete Geometrical Visualization Using the Minimum Spanning Forest of Heuristic Target Net Offset and the Topology Escape Routing Algorithm Version Four](./documentation/fig/result.png)

## Project Contents

```txt
tadpole
├── CMakeLists.txt
├── Io_drc
├── README.html
├── README.md
├── bin
├── build
├── documentation
├── inc
├── result
├── src
├── tags
└── via_info
```

- 首先我們看到 CMakeLists.txt 檔案，此 project 以 cmake 編譯，其依賴的 library 有:
  - [boost](https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/index.html): 我會需要裡面的 boost graph library, geometry library polygon library
  - [fmt](https://github.com/fmtlib/fmt): 有些 cout 是用這個 library 印出的
  - [xlsxwriter](http://libxlsxwriter.github.io): 忘記哪裡用到了，說不定不用
  - [gnuplot 5.2](http://www.gnuplot.info): 我的 geometrical visualization 會輸出一個 gnuplot 5 語法的 script，要更新到 gnuplot 5  
  - [python3 library openpyxl](https://openpyxl.readthedocs.io/en/stable/):
    這是用在最後將所有結果自動化產生一個 xlsx 檔，./src/experimentw.cpp 及 ./src/write_to_xlsx.py 會使用到，產生的結果請看 result 資料夾的內容。

    建議細看我的 cmake 的寫法，大概的重點是：

  - cpp 檔都放在 src 資料夾
  - header 檔都放在 inc 資料夾
  - tags 檔案是由 ctags 產生，trace code 會需要，強烈建議要用。
  - Input 資料都放在 Io_drc 及 via_info 資料夾裡
  - 編譯時請將 pwd 設定在 tadpole，然後使用以下 command 編譯：

``` bash
cmake . -B build
make -C build
```

- cmake 自動產生的文件會放在 build 資料夾裡，產生的執行檔會放在 bin 資料夾裡，請將 pwd 設為 tadpole，然後直接執行執行檔就可以了。e.g.

``` bash
./bin/demo_1
./bin/experiment2
```

- 我的實作內容都放在 ./inc/topology.h 檔案裡。所有用到我的實作的 source file  都必須包含 `#include <topology.h>`。其具體如何使用，請看 ./src/experiment2.cpp，這個 source file 就是我的論文中整個 topology escape routing flow 的實作，要 trace code 可以從這裡開始。

## Developer Notes

- 若要了解我的論文中 triple list table 怎麼使用，可以見 demo_1.cpp ～ demo_6.cpp，這幾個 demo 是用來產生我論文中的示意圖的，可以自己玩一下，有助於了解 triple list 的核心實作。若遇到 segmentation fault 或是很奇怪的結果，就是你 make_slice 了不在相同 slice 的 topology vertices，這樣是不合法的所以產生未定義行為。唯注意 demo 中 `syc::topology::model::v_01::sketchable_forest<>` 這個 class 就是 triple list table 的實作，v_01 代表是第一版，名字不一樣是因為後來論文的最終名字與一開始取的名字不一樣。

- 我們來看看 ./inc/topology.h 。建議使用 ctags 來 trace code，可以看到右邊的 panel，這是 [vim-ctags plugin](https://github.com/webastien/vim-ctags) 顯示 ctags  產生的 tags 檔的結果。![Screen shot 1](./documentation/fig/Screen_shot_1.png) 我的實作都放在 syc::topology  namespace 之下，v_01 代表第一版，sketchable_forest就是實作在這一版本裡，第二第三版都會繼承這個 class，而第二版是過渡期，不用進去看，需要看 v_01 及 v_03 的內容就可以了

- 其他沒提到的檔案都是過度檔案，有些 bug 不用理它們，但也可以進去看看參考一下。因為趕畢業，很多寫得很亂抱歉了ＱＱ
