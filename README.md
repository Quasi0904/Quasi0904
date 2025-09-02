代码用于处理HST的图像，对M31的星团进行测光并绘制CMD图像及等年龄线。

主要步骤如下：

⭐根据Johnson L.C.的星表 PHAT Stellar Cluster Survey II. Andromeda Project Cluster Catalog找到需要研究的星团的相关信息，主要是APID,RADeg,DECDeg,Rap。
使用hst的网站https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html找到包含该星团的图像并下载。下载的波段是F814W（对应文件夹jbf*10）、F475W（对应文件夹jbf*10），储存位置是'.../m31/origion/HST/jbf*'。

⭐运行preprocess.py代码：
    出现“请输入原始数据文件夹名前缀 (例如: jbf8010)：”，输入星团对应图像的文件夹名前缀；
    如果出现“警告: 输出目录已存在且不为空”，则包含该星团的图像已经经过预处理，可以跳过后面的步骤；
    输入后，preprocess.py会调用DOLPHOT对图像进行预处理：acsmask、splitgroups、calcsky;
    处理完的图像储存在'.../m31/origion/raw/jbf*'

⭐运行photometry.py代码：
    出现“请输入原始数据文件夹名除最后两位10或20的前缀：”，输入星团对应图像的文件夹名前缀；
    出现“请输入文件夹ID (例如: 1 对应 J1): ”，输入星团的APID；
    出现“请输入目标赤经 (deg): ”，输入星团的RA；
    出现“请输入目标赤纬 (deg):”，输入星团的DEC；
    出现“请输入有效半径 (arcsec)”，输入星团的Rap；
    出现“配置信息确认:”，检查配置信息是否正确，储存位置是否正确，如确定无误，输入“y”。
    输入和保存目录为
      原始数据目录: '.../m31/origion/raw/jbf*'
      输出目录: '.../m31/JHONSON/J{APID}/814(475)'
    这个代码会根据输入的星团坐标和有效半径剪裁HST图像，并自动生成PSF测光所需参数的.param文件，调用DOLPHOT进行测光，并提取DOLPHOT输出星表中有效源（不过曝、判定为恒星），并把这些源的信息转为jh2000格式，储存在'.../m31/JHONSON/J{APID}/J{APID}_814(475).txt'中

⭐使用TOPCAT匹配photometry.py输出的J"APID"_814.txt和J"APID"_475.txt文件中的源，提取出在两个波段中都被识别出并测光的源。保存为"APID".txt,保存在 '.../m31/JHONSON/J"APID"'

⭐运行CMD.py代码：
    根据研究的星团，修改代码中“APID”，代码会自动绘制该星团的CMD图像，可以根据需要自行修改上下限。
    图像会被保存在.../m31/JHONSON/J{APID}/{APID}.png

