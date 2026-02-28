#from PySide6 import QtCore, QtWidgets, QtGui
#from PyQt6 import  QtCore
from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QWidget

import sys


class AnotherWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        self.label = QLabel("Another Window")
        layout.addWidget(self.label)
        self.setLayout(layout)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.button = QPushButton("Push for Window") 
        self.text = QLabel("Window 2")
        self.setCentralWidget(self.button)
        #self.setWindowTitle("My first app")
        self.resize(800,600)
        #self.layout = QVBoxLayout(self)
        #self.layout.addWidget(self.button)
        #self.layout.addWidget(self.text)

        self.button.clicked.connect(self.show_new_window)
        #self.setLayout(self.layout)

    def show_new_window(self,checked):
        w = AnotherWindow()
        w.show()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()