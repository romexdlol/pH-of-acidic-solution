#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 20:13:53 2026

@author: terrylin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 20:13:53 2026

@author: terrylin
"""

import numpy as np


class AcidDatabase(object):
    
    def __init__(self, file_name = 'acids.txt'):
        '''
        

        Parameters
        ----------
        file_name : string
            DESCRIPTION. The default is 'acids.txt'.

        Returns
        -------
        None.

        '''
        
        self.filename = file_name
        self.acids = self.load_acids(file_name)
  
    def load_acids(self, file_name):
        '''
        

        Parameters
        ----------
        file_name : string
            DESCRIPTION.

        Returns
        -------
        acids : dictionary
            Maps acid chemical formulas to a list of their pka values.

        '''
        acids = {}
        file = open(file_name,"r")
        
        for line in file:
            line = line.strip()
            if line != '':
                parts = line.split(':')
                formula = parts[0].strip()
                pkas = parts[1].split(',')
                floated_pkas = []
                for pka in pkas:
                    floated_pkas.append(float(pka))
                acids[formula] = floated_pkas
        
        file.close()
        return acids
    def has_acid (self, c_form):
        '''
        

        Parameters
        ----------
        c_form : string
            Chemical formula of the acid.

        Returns
        -------
        Boolean
            True if the chemical formula of an acid is in the acid database. False otherwise

        '''
        return c_form in self.acids
    
    def get_acids(self):
        '''
        

        Returns
        -------
        keys : list
            List of the chemical formulas of the acids in the database.

        '''
        keys = self.acids.keys()
        return keys
    
    def get_pka (self, c_form):
        '''
        

        Parameters
        ----------
        c_form : string
            Chemical formula of the acid.

        Raises
        ------
        ValueError
            Error occurs if the acid is not in the database.

        Returns
        -------
        list
            List of pka values of the acid with chemcial formula of c_form.

        '''
        if not self.has_acid(c_form):
            raise ValueError('Acid not Found')
        return self.acids[c_form]
    
    def add_acid (self, c_form, pka_list):
        '''
        

        Parameters
        ----------
        c_form : string
            Chemical formula of the acid.
        pka_list : list
            list of pka values associated with an acid with chemcial formula c_form.

        Returns
        -------
        None.

        '''
        self.acids[c_form] = pka_list
        pka_stringed = []
        for pka in pka_list:
            pka_stringed.append(str(pka))  
    
        pka_line = ','.join(pka_stringed)   
    
        with open(self.filename, "a") as f:
            f.write(c_form + ':' + pka_line + '\n')   
            
    

class AcidSolution (object):
    
    def __init__(self, c_formula, conc, pka_list):
        '''
        

        Parameters
        ----------
        c_formula : string
            Chemical formula of an acid.
        conc : float
            Molar concentration of acid solution.
        pka_list : list
            List of pka values associated with a specific acid.

        Returns
        -------
        None.

        '''
        self.c_formula = c_formula
        self.conc = conc
        self.pka_list = pka_list
        
    def get_proticity (self):
        '''
        

        Returns
        -------
        int
            Returns the proticity (number of acidic protons) in the acid.

        '''
        return len(self.pka_list)
    
    def get_conc (self):
        '''
        

        Returns
        -------
        float
            Returns the molar concentration of the acid solution.

        '''
        return self.conc
    
    def get_pka (self):
        '''
        

        Returns
        -------
        list
            Returns a list of the pka values of the acid.

        '''
        return self.pka_list
    
    def get_most_acidic (self):
        '''
        

        Returns
        -------
        float
            Returns the lowest pka value of the acid.

        '''
        return self.pka_list[0]
    
    
    
def calc_pH_strong(acidSolution):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution

    Returns
    -------
    float
        Returns the pH of the solution.

    '''
    H = acidSolution.get_conc()
    if H < 10**-6:
        return numerical(acidSolution)
    return -np.log10(H)
        
def calc_pH_weak_Monorotic(acidSolution):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution

    Returns
    -------
    float
        Returns the pH of the solution.

    '''
    
    conc = acidSolution.get_conc()
    Ka = 10** (-acidSolution.get_pka()[0])
    h_test = (Ka *conc)**0.5
    if (h_test/conc * 100 <= 5):
        if h_test > 10**-6: 
            return -np.log10(h_test)
        
    coefficients = [1, Ka, -Ka*conc]
    roots = np.roots(coefficients)
    for root in roots:
        if np.isreal(root):
            test = root.real
            if test < conc and test > 0 and test > 10**-6:
               return -np.log10(test)
    return numerical(acidSolution)


def chooseModel(acidSolution):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution

    Returns
    -------
    string
        Returns the model in which we should use to calculate the pH of the acid solution

    '''
    proticity = acidSolution.get_proticity()
    pka = acidSolution.get_pka()
    if proticity == 1:
        if pka[0] <= 0:
            return 'strong'
        else:
            return 'weak_monoprotic'
    else:
        return 'numerical'
    
def calc_pH(acidSolution):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution

    Returns
    -------
    float
        Returns the pH of the solution.

    '''
    model = chooseModel(acidSolution)
    
    if model == 'strong':
        return calc_pH_strong(acidSolution)
    elif model == 'weak_monoprotic':
        return calc_pH_weak_Monorotic(acidSolution)
    else:
        return numerical(acidSolution)
    
def alpha_fraction(acidSolution, H):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution
    H : float
        Molar concentration of H+.

    Returns
    -------
    alpha : list
        list of the fractions of specific forms of the acid at a given pH.

    '''
    pka_list = acidSolution.get_pka()
    ka_list = []
    for pka in pka_list:
        ka_list.append(10**-pka)
    ka_prod = [1.0]
    proticity = acidSolution.get_proticity()
    for i in range (proticity):
        ka_prod.append(ka_prod[-1] * ka_list[i])
    denom = 0.0
    for i in range (proticity + 1):
        denom += ka_prod[i] * H**(proticity-i)
    alpha = []
    for i in range (proticity+1):
        alpha_i = (ka_prod[i] * (H**(proticity-i)))/denom
        alpha.append(alpha_i)
    return alpha
    
def charge_balance(acidSolution, H):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution.
    H : float
        Molar concentration of H+.

    Returns
    -------
    float
        The charge imabalance at a specific pH.

    '''
    alpha = alpha_fraction(acidSolution,H)
    neg = 0.0
    conc = acidSolution.get_conc()
    for charge in range (1, len(alpha)):
        neg += charge * conc * alpha[charge]
    OH = 10**-14/H
    return H - (OH + neg)

def numerical(acidSolution):
    '''
    

    Parameters
    ----------
    acidSolution : AcidSolution.

    Raises
    ------
    ValueError
        Value Error if the numerical solver does not converge after a specific number of iterations.

    Returns
    -------
    float
        Returns the pH of the solution.

    '''
    highest = acidSolution.get_conc() * acidSolution.get_proticity()+10**-7
    lowest = 10**-7
    max_iter = 1000
    for i in range (max_iter):
        mid = (lowest + highest)/2
        guess = charge_balance(acidSolution, mid)
        if abs(guess) <= 10**-8:
            return -np.log10(mid)
        if guess >0:
            highest = mid
        else:
            lowest = mid
    raise ValueError('Did not converge')

def string_to_list_pka(string_of_pka):
    '''
    

    Parameters
    ----------
    string_of_pka : string
        string of pKa values separated by commas.

    Returns
    -------
    floated_pka : list
        list of pka values.

    '''
    listed_pka = string_of_pka.split(',')
    floated_pka=[]
    for pka in listed_pka:
        floated_pka.append(float(pka))
    return floated_pka

def test():
    db = AcidDatabase("acids.txt")
    print("pH Calculator")
    print("-------------")
    stay = 'Y'
    while stay.upper().strip() == 'Y':
        c_formula = input('What is the chemical formula of your acid: ').strip()
        conc = float(input('What is the Molar concentration of your acid: ').strip())
        if db.has_acid(c_formula):
            pka_list = db.get_pka(c_formula)
        else:
            string_of_pka = input('The database does not include your acid. Please enter the pka values of your acid from first to last pka separated by commas: ').strip()
            pka_list = string_to_list_pka(string_of_pka)
            db.add_acid(c_formula, pka_list)
        acidSolution = AcidSolution(c_formula, conc, pka_list)
        print('pH of your solution:'  + str(round(calc_pH(acidSolution), 2)))
        stay = input('Would you like to calculate the pH of another acid? (Y/N): ').strip()
        

if __name__ == "__main__":
    test()
