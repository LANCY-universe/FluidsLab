from PyQt5.QtWidgets import QMainWindow
from fluidslabui import Ui_MainWindow

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QProgressBar, QDialog, QGridLayout
from PyQt5.QtCore import QThread, pyqtSignal
import os
import random
import pandas as pd
import numpy as np
from numpy import log10,abs
from scipy.optimize import root,fsolve,leastsq
import DEW_HL
import DEW_EQ
import time

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure

class NewThread_pH(QThread):

    finishSignal = pyqtSignal(str)
    finishSignal_2 = pyqtSignal(str)
    progressBarSignal = pyqtSignal(int)

    def __init__(self,  ele_list, P, T, logfO2,set_pH,parent=None):
        super(NewThread_pH, self).__init__(parent)
        
        self.ele_list = ele_list 
        self.P = P
        self.T = T
        self.logfO2 = logfO2
        self.set_pH = set_pH
        
    def run(self):
        time_start = time.time()
        # 2、获得元素对应的所有pair_species和main_species
        pair_species_database = pd.read_excel(io='./data_species0.3.xlsx',sheet_name='pair_species',header=0,index_col='element')
        main_species_database = pd.read_excel(io='./data_species0.3.xlsx',sheet_name='main_species',header=0,index_col='element')
        pair_species_list = []
        for i in range(len(pair_species_database.columns)):
            tip = True 
            for ii in range(len(pair_species_database.iloc[:,i])):
                if self.ele_list[ii] == 0 and ii != 0:
                    if pair_species_database.iloc[:,i][ii] == 0:
                        continue
                    else:
                        tip = False  #删除这个species，因为包含了多余元素
            if tip:
                pair_species_list.append(pair_species_database.columns[i])
        main_species_list = []
        for i in range(len(main_species_database.columns)):
            tip = True 
            for ii in range(len(main_species_database.iloc[:,i])):
                if self.ele_list[ii] == 0 and ii != 0:
                    if main_species_database.iloc[:,i][ii] == 0:
                        continue
                    else:
                        tip = False  #删除这个species，因为包含了多余元素
            if tip:
                if str(main_species_database.columns[i]) == '<H+>' or str(main_species_database.columns[i]) == '<H2O>':
                    continue
                else:
                    print(type(main_species_database.columns[i]))
                    main_species_list.append(main_species_database.columns[i])
        print(main_species_list)
        # 3、获得这些pair_species对应温压下的logk
        pair_species_logk_list = []

        for i in range(len(pair_species_list)):
            pair_species_logk_list.append(DEW_EQ.logk(self.P,self.T,pair_species_list[i]))

        self.progressBarSignal.emit(0)
        # 4、构建方程
        name_unknows = main_species_list + pair_species_list
        print(name_unknows)
        # print(len(pair_species_list))
        list_file_data = open('DEW_dic_data0.3.dat','r').readlines()

        # 4.1、构建logk平衡方程
        list_species_pair_reaction_line = []
        for i in range(len(pair_species_list)):
            for ii in range(len(list_file_data)):
                if pair_species_list[i] in list_file_data[ii]:
                    list_species_pair_reaction_line.append(list_file_data[ii+5])
        list_species_pair_equations = []
        for i in range(len(list_species_pair_reaction_line)):
            line_x = list_species_pair_reaction_line[i].strip().split(";")
            # print(line_x)
            x = [ii.strip() for ii in line_x if ii !='']
            a_equation = ''
            for ii in range(len(x)):
                if x[ii].split()[0] == 'O2(G)':
                    a_equation+= '{}*logfO2'.format(x[ii].split()[1])
                if x[ii].split()[0] == 'H+':
                    a_equation+= '{}*(-set_pH)'.format(x[ii].split()[1])
                if x[ii].split()[0] == 'H2O':
                    a_equation+= '{}*0'.format(x[ii].split()[1])

                for iii in range(len(name_unknows)):
                    if '<'+x[ii].split()[0]+'>' == name_unknows[iii]:
                        a = '{}*log10(x[{}]){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        a_equation+=a
                        break

            a_equation+='-({})'.format(round(pair_species_logk_list[i],4))
            list_species_pair_equations.append((a_equation))

        # 4.2 构建质量守恒方程
        # 不包含 O
        list_mass_balance_equations = []
        for i in range(len(self.ele_list)):
            a_equation = ''
            if i == 0 or i == 7 or self.ele_list[i] == 0:
                continue
            for ii in range(len(main_species_list)):
                for iii in range(len(main_species_database[main_species_list[ii]])):
                    if main_species_database[main_species_list[ii]][i] == 0:
                        continue
                    else:
                        a_equation += '+({}*x[{}])'.format(main_species_database[main_species_list[ii]][i],ii)
                        break
            for ii in range(len(pair_species_list)):
                for iii in range(len(pair_species_database[pair_species_list[ii]])):
                    if pair_species_database[pair_species_list[ii]][i] == 0:
                        continue
                    else:
                        a_equation += '+({}*x[{}])'.format(pair_species_database[pair_species_list[ii]][i],
                                                           ii+len(main_species_list))
                        break

            a_equation += '-({})'.format(self.ele_list[i])
            list_mass_balance_equations.append(a_equation)

        # 4.3 构建电荷守恒方程组
            ## 1 charge balance equations
        list_species_charge_line = []
        for i in range(len(name_unknows)):
            for ii in range(len(list_file_data)):
                if name_unknows[i] in list_file_data[ii]:
                    list_species_charge_line.append(list_file_data[ii+1].split(';')[0].strip().split('charge')[1].strip())
        charge_equation = ''
        for i in range(len(list_species_charge_line)):
            if float(list_species_charge_line[i]) != 0:
                charge_equation += '{}*x[{}]'.format(list_species_charge_line[i],i)
        charge_equation += '+{}'.format(0.1**(set_pH))

        list_all_equations = []
        list_all_equations.extend(list_species_pair_equations)
        list_all_equations.extend(list_mass_balance_equations)
        list_all_equations.append(charge_equation)
        self.progressBarSignal.emit(1)


        # 5、解方程

        Z_list =[]
        for i in range(len(list_species_charge_line)):
            Z_list.append(float(list_species_charge_line[i]))
        print(Z_list)
        result_list = []

        print(len(list_all_equations))
        print(len(name_unknows))
        

        x0 = [0.01]*len(name_unknows)
        print(list_all_equations)
        global logr_list
        logr_list = DEW_HL.all_Debye_huckel(self.P*1000,self.T,x0,Z_list,DEW_EQ.gamma_a(self.P,self.T),DEW_EQ.gamma_b(self.P,self.T))
        
        print(logr_list)

        def f(x):
            eqs = []
            for k in range(len(list_all_equations)) :
                eval('eqs.append({})'.format(list_all_equations[k]))
            return eqs          
        x0 = list(leastsq(f,x0)[0])
        result_list.append(x0)
      
        for i in range(6):
            logr_list = DEW_HL.all_Debye_huckel(self.P*1000,self.T,x0,Z_list,DEW_EQ.gamma_a(self.P,self.T),DEW_EQ.gamma_b(self.P,self.T))
            x0 = list(leastsq(f,x0)[0])
            result_list.append(x0)
            self.progressBarSignal.emit(i+2)
        print(result_list[-1])
        
        fw = open('./0_species0.3.txt','w')
        screen_str_species = 'Name\tmole/L\tcharge\n'
        for i in range(len(result_list[-1])):
            fw.write(name_unknows[i]+'   '+str(result_list[-1][i])+ '    '+ str(list_species_charge_line[i])+'\n')
            screen_str_species+=(name_unknows[i]+'\t'+str(round(result_list[-1][i],6))+ '\t'+ str(list_species_charge_line[i])+'\n')
        fw.write(str(self.P)+'\n')
        fw.write(str(self.T)+'\n')
        fw.close()
        time_end = time.time()



        # 检查矿物沉淀函数

        # 1、读输入文件，元素，温压
        # 2、读species含量
        species_name_list = name_unknows
        species_concentration_list = result_list[-1]



        # 3、判断会出现的minerals
        minerals_database = pd.read_excel(io='./data_species0.3.xlsx',sheet_name='minerals',header=0,index_col='element')
        minerals_list = []
        for i in range(len(minerals_database.columns)):
            tip = True 
            for ii in range(len(minerals_database.iloc[:,i])):
                if self.ele_list[ii] == 0:
                    if minerals_database.iloc[:,i][ii] == 0:
                        continue
                    else:
                        tip = False  #删除这个minerals，因为包含了多余元素
            if tip:
                minerals_list.append(minerals_database.columns[i])
        print(minerals_list)

        # 4、找到这些minerals对应温压的logk
        minerals_logk_list = []

        for i in range(len(minerals_list)):
            minerals_logk_list.append(DEW_EQ.logQ_mineral(self.P,self.T,minerals_list[i]))
        print(minerals_logk_list)
        # 5、构造计算矿物logk的公式
        list_file_data = open('DEW_dic_data0.3.dat','r').readlines()
        minerals_reaction_line = []
        for i in range(len(minerals_list)):
            for ii in range(len(list_file_data)):
                if minerals_list[i] in list_file_data[ii]:
                    minerals_reaction_line.append(list_file_data[ii+4])
                    print(minerals_reaction_line)
        list_minerals_equations = []
        for i in range(len(minerals_reaction_line)):
            line_x = minerals_reaction_line[i].strip().split(";")
            print(line_x)
            x = [ii.strip() for ii in line_x if ii !='']

            a_equation = ''
            for ii in range(len(x)):
                # if x[ii].split()[0] == 'H2O':
                #   a_equation+= '{}*log10(aH2O)'.format(x[ii].split()[1])
                if x[ii].split()[0] == 'O2(G)':
                    a_equation+= '{}*logfO2'.format(x[ii].split()[1])                
                for iii in range(len(species_name_list)):
                    if '<'+x[ii].split()[0]+'>' == species_name_list[iii]:
                        if x[ii].split()[0] == 'H2O':
                            a = '{}*log10(x[{}]/55.55){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        else:
                            a = '{}*log10(x[{}]){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        a_equation+=a
                        break
            a_equation+='-({})'.format(round(minerals_logk_list[i],4))
            list_minerals_equations.append((a_equation))
        print(list_minerals_equations)
        ##################################3
        ############  删除Gamma
        ###################################
        
        # # 7、读取gammaAB
        # gamma_sheet = pd.read_excel('gamma.xlsx',sheet_name='gamma')
        # P_list = gamma_sheet['P/kbar']
        # T_list = gamma_sheet['T/C']
        # a_gamma = gamma_sheet['a_gamma']
        # b_gamma = gamma_sheet['b_gamma']
        self.progressBarSignal.emit(7)
        # 算affinity
        def affinity(x,list_minerals_equations_aff):
            # print('aaaa')
            # print(x)
            '''
            x是species的浓度列表,list_minerals_equations_aff是矿物反应
            '''

            logr_list = DEW_HL.all_Debye_huckel(self.P*1000,self.T,x,Z_list,DEW_EQ.gamma_a(self.P,self.T),DEW_EQ.gamma_b(self.P,self.T))
            affinity_list = []
            for i in range(len(list_minerals_equations_aff)):
                affinity_list.append(eval('({})'.format(list_minerals_equations[i])))
            print(affinity_list)
            return affinity_list
        
        self.progressBarSignal.emit(8)
        # 输出 是否饱和str
        affinity_minerals_list = affinity(species_concentration_list,list_minerals_equations)
        saturated_mineral_str = ''
        for i in range(len(affinity_minerals_list)):
            if affinity_minerals_list[i] >= 0:
                saturated_mineral_str+=('* ' +minerals_list[i]+' (Saturated) '+'\n')
            else:
                saturated_mineral_str+=( minerals_list[i]+' (Unsaturated) '+'\n')
        self.progressBarSignal.emit(0)
        self.finishSignal.emit(screen_str_species)
        self.finishSignal_2.emit(saturated_mineral_str)
        print(time_end - time_start)
        return


class NewThread(QThread):

    finishSignal = pyqtSignal(str)
    finishSignal_2 = pyqtSignal(str)
    progressBarSignal = pyqtSignal(int)

    def __init__(self,  ele_list, P, T, logfO2,parent=None):
        super(NewThread, self).__init__(parent)
        
        self.ele_list = ele_list 
        self.P = P
        self.T = T
        self.logfO2 = logfO2
        
    def run(self):
        time_start = time.time()
        # 2、获得元素对应的所有pair_species和main_species
        pair_species_database = pd.read_excel(io='./data_species0.3.xlsx',sheet_name='pair_species',header=0,index_col='element')
        main_species_database = pd.read_excel(io='./data_species0.3.xlsx',sheet_name='main_species',header=0,index_col='element')
        pair_species_list = []
        for i in range(len(pair_species_database.columns)):
            tip = True 
            for ii in range(len(pair_species_database.iloc[:,i])):
                if self.ele_list[ii] == 0 and ii != 0:
                    if pair_species_database.iloc[:,i][ii] == 0:
                        continue
                    else:
                        tip = False  #删除这个species，因为包含了多余元素
            if tip:
                pair_species_list.append(pair_species_database.columns[i])
        main_species_list = []
        for i in range(len(main_species_database.columns)):
            tip = True 
            for ii in range(len(main_species_database.iloc[:,i])):
                if self.ele_list[ii] == 0 and ii != 0:
                    if main_species_database.iloc[:,i][ii] == 0:
                        continue
                    else:
                        tip = False  #删除这个species，因为包含了多余元素
            if tip:
                main_species_list.append(main_species_database.columns[i])
        # 3、获得这些pair_species对应温压下的logk
        pair_species_logk_list = []

        for i in range(len(pair_species_list)):
            pair_species_logk_list.append(DEW_EQ.logk(self.P,self.T,pair_species_list[i]))
        self.progressBarSignal.emit(0)
        # 4、构建方程
        name_unknows = main_species_list + pair_species_list
        print(name_unknows)
        # print(len(pair_species_list))
        list_file_data = open('DEW_dic_data0.3.dat','r').readlines()

        # 4.1、构建logk平衡方程
        list_species_pair_reaction_line = []
        for i in range(len(pair_species_list)):
            for ii in range(len(list_file_data)):
                if pair_species_list[i] in list_file_data[ii]:
                    list_species_pair_reaction_line.append(list_file_data[ii+5])
        list_species_pair_equations = []
        for i in range(len(list_species_pair_reaction_line)):
            line_x = list_species_pair_reaction_line[i].strip().split(";")
            # print(line_x)
            x = [ii.strip() for ii in line_x if ii !='']
            a_equation = ''
            for ii in range(len(x)):
                if x[ii].split()[0] == 'O2(G)':
                    a_equation+= '{}*logfO2'.format(x[ii].split()[1])
                for iii in range(len(name_unknows)):
                    if '<'+x[ii].split()[0]+'>' == name_unknows[iii]:
                        if x[ii].split()[0] == 'H2O':
                            a = '{}*log10(x[{}]/55.55){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        else:
                            a = '{}*log10(x[{}]){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        a_equation+=a
                        break
            a_equation+='-({})'.format(round(pair_species_logk_list[i],4))
            list_species_pair_equations.append((a_equation))

        # 4.2 构建质量守恒方程
        # 不包含 O
        list_mass_balance_equations = []
        for i in range(len(self.ele_list)):
            a_equation = ''
            if i != 0:
                if i == 7 or self.ele_list[i] == 0:
                    continue
            for ii in range(len(main_species_list)):
                for iii in range(len(main_species_database[main_species_list[ii]])):
                    if main_species_database[main_species_list[ii]][i] == 0:
                        continue
                    else:
                        a_equation += '+({}*x[{}])'.format(main_species_database[main_species_list[ii]][i],ii)
                        break
            for ii in range(len(pair_species_list)):
                for iii in range(len(pair_species_database[pair_species_list[ii]])):
                    if pair_species_database[pair_species_list[ii]][i] == 0:
                        continue
                    else:
                        a_equation += '+({}*x[{}])'.format(pair_species_database[pair_species_list[ii]][i],
                                                           ii+len(main_species_list))
                        break
            if i == 0:
                a_equation += '-({})-111.111'.format(self.ele_list[i])
            else:
                a_equation += '-({})'.format(self.ele_list[i])
            list_mass_balance_equations.append(a_equation)

        # 4.3 构建电荷守恒方程组
            ## 1 charge balance equations
        list_species_charge_line = []
        for i in range(len(name_unknows)):
            for ii in range(len(list_file_data)):
                if name_unknows[i] in list_file_data[ii]:
                    list_species_charge_line.append(list_file_data[ii+1].split(';')[0].strip().split('charge')[1].strip())
        charge_equation = ''
        for i in range(len(list_species_charge_line)):
            if float(list_species_charge_line[i]) != 0:
                charge_equation += '{}*x[{}]'.format(list_species_charge_line[i],i)

        list_all_equations = []
        list_all_equations.extend(list_species_pair_equations)
        list_all_equations.extend(list_mass_balance_equations)
        list_all_equations.append(charge_equation)
        self.progressBarSignal.emit(1)


        # 5、解方程

        Z_list =[]
        for i in range(len(list_species_charge_line)):
            Z_list.append(float(list_species_charge_line[i]))

        result_list = []

        x0 = [0.01]*len(list_all_equations)
        x0[name_unknows.index('<H2O>')] = 55.55
        global logr_list
        logr_list = DEW_HL.all_Debye_huckel(self.P*1000,self.T,x0,Z_list,DEW_EQ.gamma_a(self.P,self.T),DEW_EQ.gamma_b(self.P,self.T))
        
        # print(logr_list)
        # print(list_all_equations)
        def f(x):
            eqs = []
            for k in range(len(list_all_equations)) :
                eval('eqs.append({})'.format(list_all_equations[k]))
            return eqs          
        x0 = list(leastsq(f,x0)[0])
        result_list.append(x0)

        for i in range(6):
            logr_list = DEW_HL.all_Debye_huckel(self.P*1000,self.T,x0,Z_list,DEW_EQ.gamma_a(self.P,self.T),DEW_EQ.gamma_b(self.P,self.T))
            x0 = list(leastsq(f,x0)[0])
            result_list.append(x0)
            self.progressBarSignal.emit(i+2)
        print(result_list[-1])
        
        fw = open('./0_species0.3.txt','w')
        screen_str_species = 'Name\tmole/L\tcharge\n'
        for i in range(len(result_list[-1])):
            fw.write(name_unknows[i]+'   '+str(result_list[-1][i])+ '    '+ str(list_species_charge_line[i])+'\n')
            screen_str_species+=(name_unknows[i]+'\t'+str(round(result_list[-1][i],6))+ '\t'+ str(list_species_charge_line[i])+'\n')
        fw.write(str(self.P)+'\n')
        fw.write(str(self.T)+'\n')
        fw.close()
        time_end = time.time()



        # 检查矿物沉淀函数

        # 1、读输入文件，元素，温压
        # 2、读species含量
        species_name_list = name_unknows
        species_concentration_list = result_list[-1]



        # 3、判断会出现的minerals
        minerals_database = pd.read_excel(io='./data_species0.3.xlsx',sheet_name='minerals',header=0,index_col='element')
        minerals_list = []
        for i in range(len(minerals_database.columns)):
            tip = True 
            for ii in range(len(minerals_database.iloc[:,i])):
                if self.ele_list[ii] == 0:
                    if minerals_database.iloc[:,i][ii] == 0:
                        continue
                    else:
                        tip = False  #删除这个minerals，因为包含了多余元素
            if tip:
                minerals_list.append(minerals_database.columns[i])
        print(minerals_list)

        # 4、找到这些minerals对应温压的logk
        minerals_logk_list = []

        for i in range(len(minerals_list)):
            minerals_logk_list.append(DEW_EQ.logQ_mineral(self.P,self.T,minerals_list[i]))
        print(minerals_logk_list)
        # 5、构造计算矿物logk的公式
        list_file_data = open('DEW_dic_data0.3.dat','r').readlines()
        minerals_reaction_line = []
        for i in range(len(minerals_list)):
            for ii in range(len(list_file_data)):
                if minerals_list[i] in list_file_data[ii]:
                    minerals_reaction_line.append(list_file_data[ii+4])
                    print(minerals_reaction_line)
        list_minerals_equations = []
        for i in range(len(minerals_reaction_line)):
            line_x = minerals_reaction_line[i].strip().split(";")
            print(line_x)
            x = [ii.strip() for ii in line_x if ii !='']

            a_equation = ''
            for ii in range(len(x)):
                # if x[ii].split()[0] == 'H2O':
                #   a_equation+= '{}*log10(aH2O)'.format(x[ii].split()[1])
                if x[ii].split()[0] == 'O2(G)':
                    a_equation+= '{}*logfO2'.format(x[ii].split()[1])                
                for iii in range(len(species_name_list)):
                    if '<'+x[ii].split()[0]+'>' == species_name_list[iii]:
                        if x[ii].split()[0] == 'H2O':
                            a = '{}*log10(x[{}]/55.55){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        else:
                            a = '{}*log10(x[{}]){}*(logr_list[{}])'.format(x[ii].split()[1],iii,x[ii].split()[1],iii)
                        a_equation+=a
                        break
            a_equation+='-({})'.format(round(minerals_logk_list[i],4))
            list_minerals_equations.append((a_equation))
        print(list_minerals_equations)


        self.progressBarSignal.emit(7)
        # 算affinity
        def affinity(x,list_minerals_equations_aff):
            # print('aaaa')
            # print(x)
            '''
            x是species的浓度列表,list_minerals_equations_aff是矿物反应
            '''

            logr_list = DEW_HL.all_Debye_huckel(self.P*1000,self.T,x,Z_list,DEW_EQ.gamma_a(self.P,self.T),DEW_EQ.gamma_b(self.P,self.T))
            affinity_list = []
            for i in range(len(list_minerals_equations_aff)):
                affinity_list.append(eval('({})'.format(list_minerals_equations[i])))
            print(affinity_list)
            return affinity_list
        
        self.progressBarSignal.emit(8)
        # 输出 是否饱和str
        affinity_minerals_list = affinity(species_concentration_list,list_minerals_equations)
        saturated_mineral_str = ''
        for i in range(len(affinity_minerals_list)):
            if affinity_minerals_list[i] >= 0:
                saturated_mineral_str+=('* ' +minerals_list[i]+' (Saturated) '+'\n')
            else:
                saturated_mineral_str+=( minerals_list[i]+' (Unsaturated) '+'\n')
        self.progressBarSignal.emit(0)
        self.finishSignal.emit(screen_str_species)
        self.finishSignal_2.emit(saturated_mineral_str)
        print(time_end - time_start)
        return


class Figure_Canvas(FigureCanvas):
    """
    主界面创建画板
    """
    def __init__(self, width=4.2, height=3.7):
        self.fig = Figure(figsize=(width, height), dpi=70)
        super(Figure_Canvas, self).__init__(self.fig)
    # 画pie图
    def plot_pie(self):
        species_name = []
        species_con = []
        fr = open('./0_species0.3.txt','r').readlines()
        for i in range(len(fr)-2):
            if fr[i].split()[0] == '<H2O>':
                continue
            species_name.append(fr[i].split()[0])
            species_con.append(fr[i].split()[1])
        self.fig.clf()  # 清除画布
        self.ax = self.fig.add_subplot(111)  # 111表示1行1列，第一张曲线图
        self.ax.pie(species_con,labels=species_name,autopct='%.2f%%',radius=1.1)
        self.ax.legend(loc=2,bbox_to_anchor=(-0.1,1.1),fontsize=12.5)  # 添加图例


class NewWindow(QDialog):
    # 主界面画饼状图
    def __init__(self):
        super().__init__()
        self.initUI()
        self.plot_in_convas()
    def initUI(self):
        self.setWindowTitle('plot results')
        self.resize(700,700)
        self.groupBox = QtWidgets.QGroupBox(self)
        self.groupBox.setGeometry(QtCore.QRect(50, 50, 600, 600))
        self.groupBox.setObjectName("groupBox")
    def plot_in_convas(self):
        self.PieFigure = Figure_Canvas()
        self.toolbar = NavigationToolbar(self.PieFigure, self)
        self.PieFigureLayout = QGridLayout(self.groupBox)
        self.PieFigureLayout.addWidget(self.PieFigure)
        self.PieFigure.plot_pie()



class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        # self.thread = NewThread()
        self.pushButton_10.clicked.connect(self.openfile)
        self.pushButton_21.clicked.connect(self.savefile)
        self.pushButton_22.clicked.connect(self.start_calculate)
        self.pushButton_26.clicked.connect(self.plot_results)


    def openfile(self):
        # print('ccc')
        filename=QFileDialog.getOpenFileName(None,'选择输入文件',os.getcwd(), "All Files(*);;Text Files(*.txt)")
        print(list(filename))
        fr = open(str(list(filename)[0]),'r').readlines()
        element_list = eval(fr[0])
        print(element_list)
        self.lineEdit.setText(str(element_list[0])) # H
        self.lineEdit_22.setText(str(element_list[1])) # He
        self.lineEdit_2.setText(str(element_list[2])) # Li
        self.lineEdit_6.setText(str(element_list[3])) # Be
        self.lineEdit_11.setText(str(element_list[4])) # B
        self.lineEdit_10.setText(str(element_list[5])) # C
        self.lineEdit_19.setText(str(element_list[6])) # N
        self.lineEdit_18.setText(str(element_list[7])) # O
        self.lineEdit_16.setText(str(element_list[8])) # F
        self.lineEdit_15.setText(str(element_list[9])) # Ne
        self.lineEdit_3.setText(str(element_list[10])) # Na
        self.lineEdit_7.setText(str(element_list[11])) # Mg
        self.lineEdit_13.setText(str(element_list[12])) # Al
        self.lineEdit_9.setText(str(element_list[13])) # Si
        self.lineEdit_21.setText(str(element_list[14])) # P
        self.lineEdit_17.setText(str(element_list[15])) # S
        self.lineEdit_14.setText(str(element_list[16])) # Cl
        self.lineEdit_20.setText(str(element_list[17])) # Ar
        self.lineEdit_4.setText(str(element_list[18])) # K
        self.lineEdit_5.setText(str(element_list[19])) # Ca
        self.lineEdit_12.setText(str(fr[1].strip())) # P
        self.lineEdit_8.setText(str(fr[2]).strip()) # T
        self.lineEdit_23.setText(str(fr[3]).strip()) # logfo2

    def savefile(self):
        directory = QtWidgets.QFileDialog.getExistingDirectory(None,"选取文件夹","C:/")  # 起始路径
        rand_int = random.randint(0,100)
        file_directory = str(directory)+'/input{}.txt'.format(str(rand_int))       
        fw = open(file_directory,'w')
        fw.write('['+self.lineEdit.text()+','+self.lineEdit_22.text()+','+self.lineEdit_2.text()+','+self.lineEdit_6.text()+','+self.lineEdit_11.text()+','
            +self.lineEdit_10.text()+','+self.lineEdit_19.text()+','+self.lineEdit_18.text()+','+self.lineEdit_16.text()+','+self.lineEdit_15.text()+','
            +self.lineEdit_3.text()+','+self.lineEdit_7.text()+','+self.lineEdit_13.text()+','+self.lineEdit_9.text()+','+self.lineEdit_21.text()+','+self.lineEdit_17.text()+','
            +self.lineEdit_14.text()+','+self.lineEdit_20.text()+','+self.lineEdit_4.text()+','+self.lineEdit_5.text()+']'+'\n')
        fw.write(self.lineEdit_12.text()+'\n')
        fw.write(self.lineEdit_8.text()+'\n')
        fw.write(self.lineEdit_23.text())
        fw.close()
        QMessageBox.information(None, "保存提示", "保存成功：input{}.txt".format(rand_int))

    def start_calculate(self):
        print('....calculating....')
        
        #print(self.lineEdit_17.text())
        # self.textBrowser_2.setText( "0.00")
        # self.textBrowser.setText( "0.00")

        # 1、读输入文件，元素，温压, 氧逸度
        ele_str_list = eval('['+self.lineEdit.text()+','+self.lineEdit_22.text()+','+self.lineEdit_2.text()+','+self.lineEdit_6.text()+','+self.lineEdit_11.text()+','
            +self.lineEdit_10.text()+','+self.lineEdit_19.text()+','+self.lineEdit_18.text()+','+self.lineEdit_16.text()+','+self.lineEdit_15.text()+','
            +self.lineEdit_3.text()+','+self.lineEdit_7.text()+','+self.lineEdit_13.text()+','+self.lineEdit_9.text()+','+self.lineEdit_21.text()+','+self.lineEdit_17.text()+','
            +self.lineEdit_14.text()+','+self.lineEdit_20.text()+','+self.lineEdit_4.text()+','+self.lineEdit_5.text()+']')
        ele_list = []
        for i in range(len(ele_str_list)):
            ele_list.append(float(ele_str_list[i]))
        P = float(self.lineEdit_12.text()) #kbar
        T = float(self.lineEdit_8.text()) #oC
        global logfO2
        logfO2 = float(self.lineEdit_23.text())
        global set_pH
        set_pH = float(self.lineEdit_24.text())        
        if self.checkBox_1.isChecked():
            print('pH format')
            self.th=NewThread_pH(ele_list,P,T,logfO2,set_pH)
        else:
            self.th=NewThread(ele_list,P,T,logfO2)

        self.th.finishSignal.connect(self.print_screen)
        self.th.finishSignal_2.connect(self.print_screen_2)
        self.th.progressBarSignal.connect(self.progressBarFunction)
        self.th.start()


    def print_screen(self,msg):
        QMessageBox.information(None, "计算提示", "计算完成！")
        self.textBrowser.setText(msg)

    def print_screen_2(self,msg):
        self.textBrowser_2.setText(msg)

    def progressBarFunction(self,num):
        self.progressBar.setValue(num)
        #self.textBrowser_2.setText( "0.00")

    def plot_results(self):
        self.chile_Win=NewWindow()
        self.chile_Win.show()
        
    def PTfO2_window(self):
        self.chile_Win_PT=NewWindow_PTPath()
        self.chile_Win_PT.show()
