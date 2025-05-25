文档编辑人：石子睿			最新修订日期：2025/05/19			联系方式：riceszr@163.com
在运行程序前，你需要做的事：
1.将所有文件（包括代码）存储在一个全英文的路径中，并确保路径中的所有文件夹不含空格（' '）

2.如果你使用的是visual studio：
请打开调试（D）->配置属性->调试->命令参数中一栏输入单引号内的内容'-k Opt-Path -x 5'
请不要遗漏输入项之间的空格
Opt-Path应为你想要使用的配置文件的路径，如果你不知道如何设置自己需要的配置文件，你可以使用默认的配置文件，它放在IGS-Data中，名称为GNSS_Option.opt
GNSS_Option.opt默认使用单北斗L1（B1I）单频伪距单点定位
示例：'-k D:\AAA-Study\SPP-PPP\LXZ-rtklib\IGS-Data\GNSS_Option.opt -x 5'，请在复制粘贴的过程中不要复制单引号（'）
接着，打开C/C++->常规->附加包含目录，在这一栏直接复制粘贴单引号内的内容'../../../src;../../../src/Nequick/lib/private;../../../src/Nequick/lib/public'
接着，打开C/C++->预处理器->预处理器定义，在这一栏直接复制粘贴以下内容：
--------------分割线---------------请勿复制--------------------------
WIN32
NDEBUG
_CONSOLE
_CRT_SECURE_NO_DEPRECATE
ENAGLO
ENACMP
ENAGAL
NFREQ=6
TRACE
--------------分割线---------------请勿复制--------------------------
如果你使用的不是是visual studio，请在网上自行搜索以上项的输入方法。

3.以下为主要文件夹的介绍：
app：存放工程文件
src：存放软件的核心源代码
data：存放部分输入数据文件，需要注明的是BRD40650.25P是RINEX4.0版本的广播星历，其余文件均为RINEX3.0
IGS-Data：存放部分输入数据文件（推荐），其中包含几个特殊文件，以下为介绍：
crx2rnx_ALL.bat		：一个批处理文件，他可以一键将你从网上下载的'.crx'文件转换为rinex格式（默认后缀'.25o'），并将原始文件移入RawFile文件夹中，需要当前目录下存在'crx2rnx.exe'	
crx2rnx.exe			：crx2rnx程序本体
AfterSol.bat		：一个批处理文件，他可以一键删除程序运行产生的'.trace','.pospsu_res'&'.pospsu_snr'文件，并将定位结果'.pos'文件移入Solution文件夹中
XYZ2BLH.m			：一个matlab程序，实现XYZ和BLH的相互转换，并将转换结果输出在Trans-result中，详情见注释
Trans-result.txt	：XYZ2BLH.m输出结果
SOL文件夹			：存放国际数据中心的结果文件，用于与程序定位结果做参考

4.opt文件设置介绍：
以下为一定需要你修改的设置项------------------------------------------
-outfile			：在此处输入你的输出文件名称和路径
-infile				：在此处输入你的观测数据文件或星历文件路径
-infile				：在此处输入你的星历文件或观测数据文件路径
两个'-infile'无顺序要求，只要确保一个为观测数据文件，一个为星历文件即可
示例：
-outfile             =D:\AAA-Study\SPP-PPP\LXZ-rtklib\IGS-Data\WUH200CHN_R_20251210000_01D_30S_MO.pos
-infile              =D:\AAA-Study\SPP-PPP\LXZ-rtklib\IGS-Data\WUH200CHN_R_20251210000_01D_30S_MO.25o
-infile              =D:\AAA-Study\SPP-PPP\LXZ-rtklib\IGS-Data\BRDC00IGS_R_20251210000_01D_MN.rnx
以下为可能需要你修改的设置项------------------------------------------
pos1-frequency		：解算频点（详情见文档尾部）
pos1-exclsats		：被排除的卫星（不参与解算），此处解释以下为什么默认模式要排除那些卫星，C01-C18是北斗二号卫星，其余的为文档编辑时（20250519）卫星健康状态不好的卫星
pos1-navsys			：卫星系统（掩位码模式），1:GPS+2:sbas+4:GLONASS+8:Galileo+16:QZSS+32:BDS，如果你想使用多系统解算，请直接将对应数字相加，例如：GPS+BDS=1+32=33
pos1-ionoopt		：电离层模型，需要根据解算频点调整，详情见文档尾部
ant1-pos1			：测站大致位置，X方向，可以不用特别精确；同时也是输出文件中ENU结果的参考值
ant1-pos2			：测站大致位置，Y方向，可以不用特别精确；同时也是输出文件中ENU结果的参考值
ant1-pos3			：测站大致位置，Z方向，可以不用特别精确；同时也是输出文件中ENU结果的参考值
ant1-antdele		：测站天线信息，东方向，在你的观测文件'ANTENNA: DELTA H/E/N'一行中可以找到，为第二项
ant1-antdeln		：测站天线信息，北方向，在你的观测文件'ANTENNA: DELTA H/E/N'一行中可以找到，为第三项
ant1-antdelu		：测站天线信息，高程方向，在你的观测文件'ANTENNA: DELTA H/E/N'一行中可以找到，为第一项
示例：
pos1-frequency       =l1
pos1-exclsats        =C01-C18 C31 C56-C58 C61-C63
pos1-navsys          =32         # (1:gps+2:sbas+4:glo+8:gal+16:qzs+32:comp)
pos1-ionoopt         =bdsk8      # (0:off,1:brdc,2:sbas,3:dual-freq,4:est-stec,5:ionex-tec,6:qzs-brdc,7:qzs-lex,8:vtec_sf,9:vtec_ef,10:gtec,11:bdsk8,12:bdssh9)
ant1-pos1            =-2267751.3025
ant1-pos2            =5009154.7915
ant1-pos3            =3221293.2180
ant1-antdele         =0
ant1-antdeln         =0
ant1-antdelu         =0.1206
以下为几个比较重要的opt设置项------------------------------------------
pos1-posmode		：定位模式，默认为单点定位'single'
pos1-elmask			：卫星高度角(deg)阈值，默认为5
pos1-tropopt		：对流层模型

5.频点解释
在解释前我们需要知道，北斗（BDS）的频点名称变换了好几次，这一点我们在查阅RINEX3.02，RINEX3.04，RINEX3.05的文档时就能发现差别
为了统一起见，此文档使用最新的RINEX3.05及之后的版本中所使用的名称，如果你使用的观测值文件为之前的版本，而你又想做特定频点的解算，请详细查看RINEX文档，特别注意
pos1列为：对应频点在opt文件中应输入的内容
--------------分割线---------------请勿复制--------------------------
频点名称			对应频率（MHz）		Rinex文件中的伪距名（Data Pilot D+P）		pos1-frequency		应当使用的电离层模型		pos1-ionoopt
B1I				1561.098			C2I			C2Q			C2X				l1					BDSK8					bdsk8
B2I（北二专属）	1207.140			C7I			C7Q			C7X				l2					XXXXX					XXXXX
B3I				1268.520			C6I			C6Q			C6X				l5					BDSK8					bdsk8
B1C				1575.420			C1D			C1P			C1X				l6					BDGIM					bdssh9
B2a				1176.450			C5D			C5P			C5X				l7					BDGIM					bdssh9
B2b				1207.140			C7D			C7P			C7Z				l8					BDGIM					bdssh9
--------------分割线---------------请勿复制--------------------------
如果你想做双频SPP，请直接将对应频点用加号链接，并将pos1-ionoopt设置为dual-freq
支持的推荐双频组合：
频点组合			pos1-frequency		pos1-ionoopt
B1I+B3I			l1+l5				dual-freq
B1C+B2a			l6+l7				dual-freq
B1C+B2b			l6+l8				dual-freq