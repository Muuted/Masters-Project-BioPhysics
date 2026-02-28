import numpy as np
from PySide6 import QtCore, QtWidgets, QtGui
import sys
import random




class MyWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.window_num = 1
        self.Window1()


    @QtCore.Slot()
    def Window1(self):
        self.hello = ["Hello World","Hallo Welt", "Hei maailma", "Hola Mundo", "Привет мир"]

        self.button = QtWidgets.QPushButton("Click me!")
        self.text = QtWidgets.QLabel("Window 1",alignment=QtCore.Qt.AlignCenter)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button)

        self.button.clicked.connect(self.magic)
    
    @QtCore.Slot()
    def Window2(self):
        self.hello = ["Hello World","Hallo Welt", "Hei maailma", "Hola Mundo", "Привет мир"]

        self.button = QtWidgets.QPushButton("Click me!")
        self.text = QtWidgets.QLabel("Window 2",alignment=QtCore.Qt.AlignCenter)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button)

        self.button.clicked.connect(self.magic)
    
    @QtCore.Slot()    
    def magic(self):
        if self.window_num == 1:
            self.Window2()
            self.window_num = 2
        
        if self.window_num == 2:
            self.Window1()
            self.window_num = 1


class MyWidget_2(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.hello = ["window 2"]

        self.button = QtWidgets.QPushButton("Click me!")
        self.text = QtWidgets.QLabel("Window 2",alignment=QtCore.Qt.AlignCenter)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button)

        self.button.clicked.connect(self.magic)

    @QtCore.Slot()    
    def magic(self,checked):
        w = Window2()
        w.show()

class Window2(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.hello = ["window 2"]

        self.button1 = QtWidgets.QPushButton("Click me!")
        self.button2 = QtWidgets.QPushButton("Click me!")
        self.text = QtWidgets.QLabel("Window 2",alignment=QtCore.Qt.AlignCenter)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)

        self.button.clicked.connect(self.magic)

    @QtCore.Slot()    
    def magic(self):
        w = MyWidget_2()
        w.show()




if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    widget = MyWidget_2()
    widget.resize(800, 600)
    widget.show()
    sys.exit(app.exec())




