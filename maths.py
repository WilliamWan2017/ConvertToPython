# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 08:03:30 2018

@author: william
"""
 
"""
                    用正则表达式实现四则运算表达式解析器
思路：
根据计算优先级，先计算内部括号里面的运算，并用计算结果的字符串形式替换原表达式，直到没有括号运算符；
然后匹配乘法和除法的运算因子，按从左往右的顺序依次更新计算结果，最后处理加减法运算。
Tips: 需要特别注意对输入的检测和修正（如多余空格、非法输入）和浮点数的匹配
"""
import re
import functools

def checkInput(formula):
    """检测输入合法与否,是否包含字母等非法字符"""
    return not re.search("[^0-9+\-*/.()\s]",formula)

def formatInput(formula):
    """标准化输入表达式，去除多余空格等"""
    formula = formula.replace(' ','')
    formula = formula.replace('++', '+')
    formula = formula.replace('+-', '-')
    formula = formula.replace('-+', '-')
    formula = formula.replace('--', '+')
    return formula

def mul_divOperation(s):
    """乘法除法运算"""
    # 1-2*-14969036.7968254
    s = formatInput(s)
    sub_str = re.search('(\d+\.?\d*[*/]-?\d+\.?\d*)', s)
    while sub_str:
        sub_str = sub_str.group()
        if sub_str.count('*'):
            l_num, r_num = sub_str.split('*')
            s = s.replace(sub_str, str(float(l_num)*float(r_num)))
        else:
            l_num, r_num = sub_str.split('/')
            s = s.replace(sub_str, str(float(l_num) / float(r_num)))
        #print(s)
        s = formatInput(s)
        sub_str = re.search('(\d+\.?\d*[*/]\d+\.?\d*)', s)

    return s


def add_minusOperation(s):
    """加法减法运算
    思路：在最前面加上+号，然后正则匹配累加
    """
    s = formatInput(s)
    s = '+' + s
    #print(s)
    tmp = re.findall('[+\-]\d+\.?\d*', s)
    s = str(functools.reduce(lambda x, y:float(x)+float(y), tmp))
    #print(tmp)
    return s

def compute(formula):
    """无括号表达式解析"""
    #ret = formula[1:-1]
    ret = formatInput(formula)
    ret = mul_divOperation(ret)
    ret = add_minusOperation(ret)
    return ret


def calc(formula):
    """计算程序入口"""
    has_parenthesise = formula.count('(')
    if checkInput(formula):
        formula = formatInput(formula)
        while has_parenthesise:
            sub_parenthesise = re.search('\([^()]*\)', formula) #匹配最内层括号
            if sub_parenthesise:
                #print(formula+"...before")
                formula = formula.replace(sub_parenthesise.group(), compute(sub_parenthesise.group()[1:-1]))
                #print(formula+'...after')
            else:
                print('没有括号了...')
                has_parenthesise = False

        ret = compute(formula)
        print('结果为：')
        return  ret

    else:
        print("输入有误！")

if __name__ == '__main__':

    #print(mul_divOperation('17.0*7/3.0*3+1/2'))
    #print(add_minusOperation('1-2-3+1-4'))
    #print(calc('(9+8)*7/(5-2)*3+1/2'))
    s = "99+(1+2=(8)-(7+9))"
    print(eval('1.5 - 2.4 * ( (60-30 +(-40/5) * (9-2*5/3 + 7 /3*99/4*2998 +10 * 568/14 )) - (-4*3)/ (16-3*2) )'))

    print(calc('1.5 - 2.4 * ( (60-30 +(-40/5) * (9-2*5/3 + 7 /3*99/4*2998 +10 * 568/14 )) - (-4*3)/ (16-3*2) )'))
    print( calc("1 - 2 * ( (60-30 +(-9-2-5-2*3-5/3-40*4/2-3/5+6*3) * (-9-2-5-2*5/3 + 7 /3*99/4*2998 +10 * 568/14 )) - (-4*3)/ (16-3*2) )"))