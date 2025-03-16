import sys, time, os, math, random, datetime, cnbiloop
from cnbiloop import BCI_tid
import numpy as np
from math import *
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QDesktopWidget
from python_client import Trigger

class visualInterface(QWidget):
    def __init__(self):
        super(visualInterface, self).__init__()
        self.initialization()

    def initialization(self):
        """Put configuration setting here"""

        "condition"
        self.nTrial = 30
        self.curTrial = 0

        "screen"
        self.screenWidth = QDesktopWidget().availableGeometry().width()
        self.screenHeight = QDesktopWidget().availableGeometry().height()
        self.screenHalf = QtCore.QPoint(self.screenWidth/2, self.screenHeight/2) 
        self.setGeometry(0, 0, self.screenWidth, self.screenHeight)
        self.setStyleSheet("background-color: black;")


        "error"
        self.errorRate = 0.3

        "positions"
        self.nRects = 10
        self.rectLength = 80
        self.rectMargin = self.rectLength/2
        self.rectangles = {
            "color_normal": QtGui.QColor(189,189,189),
            "color_target": QtGui.QColor(222,45,38),
            "position": [self.screenHalf.x() - 5*self.rectLength - 4.5*self.rectMargin, 
                        self.screenHalf.x() - 4*self.rectLength - 3.5*self.rectMargin, 
                        self.screenHalf.x() - 3*self.rectLength - 2.5*self.rectMargin,
                        self.screenHalf.x() - 2*self.rectLength - 1.5*self.rectMargin,
                        self.screenHalf.x() - 1*self.rectLength - 0.5*self.rectMargin,
                        self.screenHalf.x() + 0*self.rectLength + 0.5*self.rectMargin,
                        self.screenHalf.x() + 1*self.rectLength + 1.5*self.rectMargin,
                        self.screenHalf.x() + 2*self.rectLength + 2.5*self.rectMargin,
                        self.screenHalf.x() + 3*self.rectLength + 3.5*self.rectMargin, 
                        self.screenHalf.x() + 4*self.rectLength + 4.5*self.rectMargin], 
        }

        self.cursor = {
            "color": QtGui.QColor(49,130,189),
            "radius": 30,   
        }

        
        "window"
        self.setWindowTitle('Visual Interface - ErrP')
        self.setFocusPolicy(QtCore.Qt.WheelFocus)
        self.setFocus()

        "time"
        self.timestamp = QtCore.QTime()
        self.timestamp.start()

        "serial port"
        self.isParallel = True
        self.trigger_value = 0
        if (self.isParallel is True):
            self.parallel = Trigger('USB2LPT')
            self.parallel.init(50)

        "cnbi loop" 
        if (sys.argv[1] == '1'): # 1 if decoding mode, 0 if calibration mode (see expLauncher)
            self.isLoop = True
        else:
            self.isLoop = False ###############need this eventually

        if (self.isLoop is True):
            self.bci = BCI_tid.BciInterface() ##########################need this
        
        self.feedback_speed = 5

        self.InitializeTrial()
        
        QtCore.QCoreApplication.processEvents()
        QtCore.QCoreApplication.flush()
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.UpdateLoop) ##########potentially need, seems like makes sure checks for input triggers every 20 msec
        self.paintCycleTime = 20  # in [ms] #############potentially need this
        timer.start(self.paintCycleTime)

        self.show()

    def InitializeTrial(self):
        self.current_target = random.choice([0, 9])
        self.current_position = random.randint(0, 9)
        self.diff = np.nan
        self.feedback_target = np.nan
        while (self.current_target == self.current_position):
            self.current_position = random.randint(0, 9)

        self.is_first = True
        self.isError = False
        self.isCorrect = False
        self.feedback_mode = False
        self.feedback_count = 0

        self.sendTrigger = True
        self.isGoal = False
        self.trialInterval = random.randint(2000, 2500)

        self.timestamp.restart()
        self.updated_time = 4000

    def UpdateLoop(self):
        if (self.isLoop is True):
            self.receiveTiD()

        self.update()

    def paintEvent(self, e):
        qp = QtGui.QPainter()
        qp.begin(self)
        self.PaintInterface(qp)
        qp.end()

    def PaintInterface(self, qp):
        """Main Drawing Function"""
        qp.setRenderHint(QtGui.QPainter.Antialiasing, True)

        "fill background"
        # qp.setBrush(QtGui.QColor(0, 0, 0))
        # qp.fillRect(0, 0, self.screenWidth, self.screenHeight, QtGui.QColor(0, 0, 0))

        "draw rectangles"
        qp.setPen(self.rectangles["color_normal"])
        for i, x_position in enumerate(self.rectangles["position"]):
            if i == self.current_target: 
                qp.fillRect(x_position, self.screenHalf.y() - self.rectLength/2, self.rectLength, self.rectLength, self.rectangles["color_target"])
            else: 
                qp.drawRect(x_position, self.screenHalf.y() - self.rectLength/2, self.rectLength, self.rectLength)
        
        elapsedTime = self.timestamp.elapsed() - self.updated_time

        "draw cursor"
        if (self.is_first == True and elapsedTime > -500):
            self.is_first = False
            self.sendTrigger = True
     
        if (self.is_first == False):

            qp.setBrush(self.cursor["color"])
            if (self.feedback_mode == False):
                center = QtCore.QPoint(self.rectangles["position"][self.current_position]+self.rectLength/2, self.screenHalf.y())
            else:
                self.updated_time = self.timestamp.elapsed()
                if (self.feedback_count == 0):
                    self.sendTrigger = True
                self.feedback_count += 1
                center = QtCore.QPoint(self.rectangles["position"][self.current_position]+self.rectLength/2+(self.feedback_count*self.feedback_speed*self.diff*-1), self.screenHalf.y())

            qp.drawEllipse(center, self.cursor["radius"], self.cursor["radius"])

        if (self.feedback_mode == True and center.x() == self.targetPosition):
            self.feedback_mode = False
            self.feedback_count = 0
            self.current_position = self.feedback_target


        "send trigger"
        self.TriggerSend()

        if (elapsedTime > self.trialInterval):
            self.sendTrigger = True
            if (self.current_position == self.current_target):
                self.InitializeTrial()

            else: 
                self.curTrial = self.curTrial + 1
                self.updated_time = self.timestamp.elapsed()
                self.trialInterval = random.randint(2000, 2500)

                if (self.current_position == 0 or self.current_position == 9):
                    temp = 1.0
                else:
                    temp = random.uniform(0, 1)

                temp_sign = np.sign(self.current_target - self.current_position)
                if (temp <= self.errorRate):
                    self.diff = -1*temp_sign
                    self.current_position = self.current_position - temp_sign # error
                    self.isError = True
                    self.isCorrect = False
                else: 
                    self.diff = temp_sign
                    self.current_position = self.current_position + temp_sign # correct
                    self.isError = False
                    self.isCorrect = True
                
                if (self.current_position == self.current_target):
                    self.isGoal = True
                else: 
                    self.isGoal = False

        if (self.curTrial > self.nTrial):
            if (self.isLoop == True):
                self.sendTiD(3)
            exit()
        
    def TriggerSend(self):
        self.trigger_value = 1*self.isCorrect + 2*self.isError + 4*self.isGoal + 8*(self.is_first == False) + 16*self.is_first + 32*self.feedback_mode
        if (self.isParallel == True and self.sendTrigger == True):
            self.parallel.signal(self.trigger_value) 
            self.sendTrigger = False

    def sendTiD(self, value):
        self.bci.id_msg_bus.SetEvent(value)
        self.bci.iDsock_bus.sendall(str.encode(self.bci.id_serializer_bus.Serialize()))

    def receiveTiD(self):
        data = None
        try:
            data = self.bci.iDsock_bus.recv(512).decode("utf-8")
            self.bci.idStreamer_bus.Append(data)
        except:
            self.nS = False
            self.dec = 0
            pass
        # deserialize ID message
        if (data):
            if (self.bci.idStreamer_bus.Has("<tobiid", "/>")):
                msg = self.bci.idStreamer_bus.Extract("<tobiid", "/>")
                self.bci.id_serializer_bus.Deserialize(msg)
                self.bci.idStreamer_bus.Clear()
                tmpmsg = int(round(float(self.bci.id_msg_bus.GetEvent())))
                

                print("Received Message: ", float(self.bci.id_msg_bus.GetEvent()))
                

                if (tmpmsg == 1):
                    print('Correct Detected')
                elif (tmpmsg == 2):
                    print('Error Detected')
                    self.feedback_mode = True
                    self.feedback_target = self.current_position - 2*self.diff
                    if (self.feedback_target == -1):
                        self.feedback_target = 0
                    elif (self.feedback_target > 9):
                        self.feedback_target = 9 
                    self.targetPosition = self.rectangles["position"][self.feedback_target] + self.rectLength/2

            elif self.bci.idStreamer_bus.Has("<tcstatus","/>"):
                MsgNum = self.bci.idStreamer_bus.Count("<tcstatus")
                for i in range(1,MsgNum-1):
                    # Extract most of these messages and trash them     
                    msg_useless = self.bci.idStreamer_bus.Extract("<tcstatus","/>")

def main():
    app = QApplication(sys.argv)
    ex = visualInterface()
    sys.exit(app.exec_())

if __name__ == '__main__': 
    main()
