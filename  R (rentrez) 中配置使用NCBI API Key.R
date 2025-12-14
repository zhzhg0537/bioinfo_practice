###在 R (rentrez) 中配置使用NCBI API Key
#方法一：永久配置（推荐，一劳永逸）
#如果你经常用，建议把它写入 R 的环境配置文件（.Renviron），这样每次启动 RStudio 自动加载，不用重复输入。
usethis::edit_r_environ()
#在打开的文件（通常是一个空白文档或已有几行配置）中，新起一行，输入：ENTREZ_KEY=你的_API_Key_粘贴在这里
#保存该文件 (Ctrl+S 或 Cmd+S)。
#重启 RStudio (Session -> Restart R)。
#重启后，你的 rentrez 所有操作就会自动默认使用这个 Key 进行加速了。
#保密： API Key 相当于你的私人通行证，不要把它发在公开的代码库（如 GitHub）里，以免被滥用导致你的账号被封禁。
#频率： 虽然有了 Key，速度上限是 10 次/秒，但在写循环代码（Loop）时，如果不需要极速，为了保险起见，有时我们还是会加一个 Sys.sleep(0.1) 稍微停顿一下，做个“礼貌”的爬虫。
#方法二：临时配置（只对当前 R 窗口有效）
#每次打开 RStudio 都要运行一次。
library(rentrez)

# 将引号里的内容替换为你刚才申请到的真实 Key
set_entrez_key("你的_API_Key_粘贴在这里")

# 测试一下，看看是不是生效了（返回不报错即可）
Sys.getenv("ENTREZ_KEY")