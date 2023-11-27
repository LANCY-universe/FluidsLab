import sys
from PyQt5.QtWidgets import QApplication,QMainWindow
from main_fluidslab import MainWindow #调用生成的py文件

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ui = MainWindow()

    ui.show()
    sys.exit(app.exec_())
