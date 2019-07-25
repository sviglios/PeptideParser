# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:25:16 2019

@author: konka
"""

import os 
import requests
import time
from selenium import webdriver
from bs4 import BeautifulSoup
from Bio import ExPASy
from Bio import SwissProt
import pandas as pd

'''
TO do:
    1. Implement search through EBI and ProteomicsDB APIs
    2. Parse arguments through command line
    3. Upload to github
    4. Tests in other PCs
'''

print('\nInitializing')

class ExtractDBs:
    """Usage: ExtractDBs(accs, out_folder), where accs is the list of Uniprot accession numbers.
        Extracts paptides from PeptideAtlas Disctinct Observed Peptides table,
        and from ProteomicsDB reference peptide table, for the given accessions.
        Because ProteomicsDB uses client based javascript rendering, the ProteomicsDB
        extraction takes place through a dummy browser, downloading
        the table and processing it. This might not work (altough by iterating 10 times over the page,
        results should be obtained), but Peptide Atlas table should
        be extracted at all times. If the output shows exception for a given accession for 
        the ProteomicsDB database, try to run the script again. Alternatively, the table could be 
        downloaded manually and renamed to
        {acc}refpep.csv, and added to the folder with the rest of the files for further 
        processing."""
        
    def __init__(self, accs, out_folder):
        
        self.accs = accs
        self.out_folder = out_folder
        

    def MakeFolder(self):
        
        if self.out_folder not in os.listdir():
            os.mkdir(self.out_folder)
            

    def GetPeptidesHTML(self):
        
        options = webdriver.ChromeOptions() 
        options.add_argument("download.default_directory=C:/Downloads")
        driver = webdriver.Chrome(options=options)
        dlpath = r"C:\Users\konka\Downloads/"
        destpath = r'C:\Users\konka\Documents\Kostas\Scripts/' + self.out_folder + '/'
        ppflag = 1
        
        for acc in self.accs:
            
            if acc + 'dopAtlas.txt' in os.listdir(self.out_folder):
                print(f'\nFile for {acc} from Peptide Atlas already exists, moving on.')
            else:
                print('Requesting pages for',acc)
                
                url = requests.get("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=479&protein_name="+acc+"&action=QUERY")
                htmltext = url.text
                
                if not url.ok:
                    print('error', url.status_code)
                
                
                soup = BeautifulSoup(htmltext, 'html.parser')
                
                doptb = soup.find_all(id='getprotein_observedlist_div')
                
                print('Parsing',acc)
                
                if len(doptb) == 1:
                    doptb = doptb[0]
                    
                    doptb_data = doptb.find_all('table')
                    
                    headers = []
                    
                    table1 = doptb_data[0]
                    rt1 = table1.tr
                    
                    ch = 0
                    while ch != 'end':
                        
                        headers.append(rt1.td.text)
                        rt1 = rt1.next_sibling
                        rt1 = rt1.next_sibling
                        if rt1 is None:
                            ch = 'end'
                    
                    table2 = doptb_data[1]
                    data = []
                    
                    rt2 = table2.tr
                    rt2 = rt2.next_sibling
                    rt2 = rt2.next_sibling
                    
                    ch = 0
                    while ch != 'end':
                        
                        rowdata = []
                        celldata = rt2.td
                        
                        ch2 = 0
                        while ch2 != 'end':
                            
                            rowdata.append(celldata.text)
                            
                            celldata = celldata.next_sibling
                            
                            if celldata == '\n':
                                celldata = celldata.next_sibling                
                            if celldata is None:
                                ch2 = 'end'
                                
                        rt2 = rt2.next_sibling
                        
                        if rt2 == '\n':
                            rt2 = rt2.next_sibling            
                        if rt2 is None:
                            ch = 'end'
                        
                        data.append(rowdata)
                            
                else:
                    print("Don't know what to do, more than one tables for", acc)
                    continue
                
                #make peptide accession numbers index, drop their column
                df_dop = pd.DataFrame(data,columns=headers)
                df_dop.index = df_dop.Accession
                df_dop = df_dop.drop('Accession', axis = 1)
                
                #if you wnat human readable format
                #fh_out = open('Files/'+acc+'dopAtlas.txt','w')
                #fh_out.write(df_dop.to_string(justify='left',index=True))
                #fh_out.close()
                
                df_dop.to_csv(self.out_folder+'/'+acc+'dopAtlas.txt')            
                    
                print(acc, len(data))
                print('Finished PeptideAtlas',acc)
            
            if acc + 'refpep.csv' in os.listdir(self.out_folder):
                print(f'\nFile for {acc} from ProteomicsDB already exists, moving on.')
            else:
                try:
                    handle = ExPASy.get_sprot_raw(acc)
                    record = SwissProt.read(handle)
                except:
                    print('Not found')
                    
                for ref in record.cross_references:
                    
                    ch = 0
                    if ref[0] == 'ProteomicsDB':
                        print('\nProteomicsDB', acc, ref[1])
                        protid = ref[1]
                        ch = 1
                    
                    if ch == 1:
                        
                        end = 0
                        for i in range(10):
                            driver.get('https://www.proteomicsdb.org/proteomicsdb/#human/proteinDetails/' + protid + '/referencePeptides')
                            button = driver.find_element_by_xpath('/html/body/div[1]/div/div[3]/article/div/table/tbody/tr[4]/td/div/div/div/div[2]/div/table/tbody/tr[2]/td/div/div[1]/div/div/div/button/span[2]')
                            print('Iteration', i+1)
                            try:
                                button.click()
                                time.sleep(3)
                                end = 1
                                break
                            except:
                                time.sleep(1)
                                continue
                        
                        if end == 0:
                            print('EXCEPTION FOR', acc, protid, '\n\nRun the script again, so the extraction can be attempted again. No post-processing will occur. If you want to activate it anyway, pass arg for postprocessing\n\n')
                            ppflag = 0
                            continue
                        
                        print('Downloaded peptides for', acc)
                        
                        if acc + 'refpep.csv' not in os.listdir(self.out_folder):
                            os.rename(dlpath+'peptides.csv', destpath + acc + 'refpep.csv')
                        else:
                            print(f'File for {acc} already exists, moving on.')
                            os.remove(r"C:\Users\konka\Downloads/" + 'peptides.csv')
                        
                        print('Finished ProteomicsDB', acc, '\n\n')
                        break
            
        time.sleep(5)
        driver.quit()
        
        return ppflag
    
    
    def UniprotFind(acc):
        
        try:
            handle = ExPASy.get_sprot_raw(acc)
            record = SwissProt.read(handle)
            print('\nUniprot data extracted for', acc)
            #sequence = record.sequence    
        except:
            return None, None, None
        
        flag = 0
        for feat in record.features:
            if feat[0] == 'SIGNAL':
                end = feat[2]
                flag = 1
            if feat[0] == 'PROPEP':
                end = feat[2]
                flag = 1
                break
        
        if flag == 0:
            return record.entry_name, 0, record.sequence_length
        
        return record.entry_name, end, record.sequence_length
    
    
    def PeptideFilter(pep):
        
        pos = len(pep) - 1
        fl = 0
        
        for a in range(len(pep)):
            if len(pep) < 8 or len(pep) > 21:
                fl = 1
            if pep[a] == 'R' or pep[a] == 'K':
                if pos != a:
                    fl = 2
            if pep[a] == 'M':
                fl = 3
            if a == pos:
                if pep[a] != 'K' and pep[a] != 'R':
                    fl = 4
        
        return fl
        
    
    def ProcessTables(self):
        
        dfs = []
        
        for acc in self.accs:
            
            name, pro_end, length = ExtractDBs.UniprotFind(acc)
            if name == None:
                print('\n', acc, 'NOT FOUND in Uniprot, skiping!\n')
                continue
            else:
                print(name, pro_end, length)
            
            df_atlas = pd.read_csv(self.out_folder + '/' + acc + 'dopAtlas.txt')
            for i in range(len(df_atlas)):
                df_atlas.loc[df_atlas.index[i],'ESS'] = df_atlas.loc[df_atlas.index[i],'ESS'].split(' ')[0]    
            df_atlas = df_atlas.sort_values('ESS', ascending = False)
            
            try:
                df_pdb = pd.read_csv(self.out_folder + '/' + acc + 'refpep.csv', sep=';')
                df_pdb = df_pdb.sort_values('ANDROMEDA_SCORE', ascending = False)
                
                df_com = pd.merge(df_atlas, df_pdb, how='inner', left_on='Sequence', right_on ='PLAIN_SEQUENCE')
                df_com = df_com.sort_values(['MISSED_CLEAVAGES','ANDROMEDA_SCORE','ESS'], ascending = [True,False,False])
                
            except:
                print('\nProteomicsDB file not available for',acc)
                df_com = df_atlas
            
            df_com = df_com.loc[~df_com.Accession.duplicated(keep='first')]
            
            dec = 0
            for i in range(len(df_com)):
                i  -= dec
                ind = df_com.index[i]
                pep = df_com.loc[ind,'Sequence']
                
                flag = ExtractDBs.PeptideFilter(pep)
                if flag in [1,2,4]:
                    df_com = df_com.drop(ind)
                    print('Dropped',pep,'flag =', flag)
                    dec += 1
                elif flag == 3:
                    if len(df_com) > 5:
                        df_com = df_com.drop(ind)
                        print('Dropped',pep,'flag = ', 3)
                        dec += 1
            
            print(f'Final number of peptides for {acc} is {len(df_com)}, picking 15 or full dataframe\n')            
            df_sel = df_com.iloc[:15,:]
            
            ind = [acc+'-'+name for i in range(len(df_sel))]
            df_sel = df_sel.set_index([ind,'Accession'])    
            dfs.append(df_sel)
                
            print('\nFinished', acc)
            
        df_final = pd.concat(dfs, sort = True)
        df_final.to_excel(self.out_folder + '/' + 'FinalPeps.xlsx',index=True)
        
        fh_out = open(self.out_folder + '/' + 'FinalPeps.txt','w')
        fh_out.write(df_final.to_string(justify='left',index=True))
        fh_out.close()
        
        return dfs
    
    
    def MouseHomolog(self, dfs):
        
        print('\nFinding mouse homologs')
        ind = 0
        new_dfs = []
        
        for acc in self.accs:
            
            try:
                handle = ExPASy.get_sprot_raw(acc)
                record = SwissProt.read(handle)
                name = record.entry_name    
            except:
                print('\nNo entry for', acc, ',continuing')
                ind += 1
                continue
            
            try:
                mname = name.split('_')[0] + '_MOUSE'
                mhandle = ExPASy.get_sprot_raw(mname)
                mrecord = SwissProt.read(mhandle)
                mseq = mrecord.sequence
                print(f'\nFound mouse homolog for {name}: {mname}')
            except:
                print(f'\nNo mouse gene entry for {acc}-{name}, continuing')
                ind += 1
                continue
            
            df = dfs[ind]
            mcol = []
            
            for row in range(len(df)):
                pepseq = df.Sequence[df.index[row]]
                print(pepseq)
                if str(pepseq) in mseq:
                    mcol.append('True')
                else:
                    mcol.append('False')
            
            df['Mouse'] = mcol
            new_dfs.append(df)
            ind += 1
            
        df_final = pd.concat(new_dfs, sort = True)
        df_final.to_excel(self.out_folder + '/' + 'MouseHomologPeptides.xlsx',index=True)


### TEST ###    
obex = ExtractDBs(['P01033','P16035','P35625','Q99727'], 'PeptidesClass')

obex.MakeFolder()

ppflag = obex.GetPeptidesHTML()

if ppflag == 1:
    dfs = obex.ProcessTables()
    obex.MouseHomolog(dfs)
    print('\n\nFinished, closing down.\n')
else:
    print('\nPost processing omitted, shutting down.\n')
### END TEST ###