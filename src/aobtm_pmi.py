# -*- coding: utf-8 -*-
"""AOBTM-PMI.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/15EAtw4w3Mo1BAHAcBL4r7EbEujZDD1T4
"""

import time
import json
import pandas as pd
import numpy as np
import re, string, unicodedata
import nltk
import inflect
from google.colab import files

nltk.download('punkt')
nltk.download('stopwords')
nltk.download('wordnet')
nltk.download('words')

from bs4 import BeautifulSoup
from nltk import word_tokenize, sent_tokenize
from nltk.corpus import stopwords
from nltk.stem import LancasterStemmer, WordNetLemmatizer

from google.colab import drive
drive.mount('/content/drive')

# Commented out IPython magic to ensure Python compatibility.
# %cd ./drive/My\ Drive/Temp/AOBTM\ DATA

df = pd.read_csv ('External-Dataset.csv')

data = df['Review']
data.head()

def strip_html(text):
    soup = BeautifulSoup(str(text), "html.parser")
    return soup.get_text()

def remove_between_square_brackets(text):
    return re.sub('\[[^]]*\]', '', text)

def denoise_text(text):
    text = strip_html(text)
    text = remove_between_square_brackets(text)
    return text

def remove_non_ascii(words):
    """Remove non-ASCII characters from list of tokenized words"""
    new_words = []
    for word in words:
        new_word = unicodedata.normalize('NFKD', word).encode('ascii', 'ignore').decode('utf-8', 'ignore')
        new_words.append(new_word)
    return new_words

def to_lowercase(words):
    """Convert all characters to lowercase from list of tokenized words"""
    new_words = []
    for word in words:
        new_word = word.lower()
        new_words.append(new_word)
    return new_words

def remove_punctuation(words):
    """Remove punctuation from list of tokenized words"""
    new_words = []
    for word in words:
        new_word = re.sub(r'[^\w\s]', '', word)
        if new_word != '':
            new_words.append(new_word)
    return new_words

def replace_numbers(words):
    """Replace all interger occurrences in list of tokenized words with textual representation"""
    p = inflect.engine()
    new_words = []
    for word in words:
        if word.isdigit():
            new_word = p.number_to_words(word, andword="", threshold=1000000)
            new_words.append(new_word)
        else:
            new_words.append(word)
    return new_words

def remove_stopwords(words):
    """Remove stop words from list of tokenized words"""
    new_words = []
    for word in words:
        if word not in stopwords.words('english'):
            new_words.append(word)
    return new_words

def stem_words(words):
    """Stem words in list of tokenized words"""
    stemmer = LancasterStemmer()
    stems = []
    for word in words:
        stem = stemmer.stem(word)
        stems.append(stem)
    return stems

def lemmatize_verbs(words):
    """Lemmatize verbs in list of tokenized words"""
    lemmatizer = WordNetLemmatizer()
    lemmas = []
    for word in words:
        lemma = lemmatizer.lemmatize(word, pos='v')
        lemmas.append(lemma)
    return lemmas

def normalize(words):
    words = remove_non_ascii(words)
    words = to_lowercase(words)
    words = remove_punctuation(words)
    words = replace_numbers(words)
    words = remove_stopwords(words)
    return words

def stem_and_lemmatize(words):
    stems = stem_words(words)
    lemmas = lemmatize_verbs(words)
    return stems, lemmas

Eng_words = set(nltk.corpus.words.words())

list_docs = []
start_time = time.time()

for i in range(data.size): #data.size
    sent = data.at[i]
    sent = denoise_text(sent)
    if isinstance(sent, str):
        List_sent = []
        List_sent = sent.split('.')
        
        for sent in List_sent:
            sent = " ".join(w for w in nltk.wordpunct_tokenize(sent) if w.lower() in Eng_words or not w.isalpha())

            words = nltk.word_tokenize(sent)
            norm_words = normalize(words)
            stems, lemmas = stem_and_lemmatize(norm_words)
            list_docs.append(lemmas)
    
    if i % 100000 == 0:
        print("Till ", i, "Success!", "in--- %s seconds ---" % (time.time() - start_time))
    
print("Terminal Execution time:") 
print("--- %s seconds ---" % (time.time() - start_time))

print(len(list_docs))

with open('list_docs.json', 'w') as f:
    json.dump(list_docs, f)

f.close()

short_texts = [x for x in list_docs if (x != [] and len(x) != 1)]
print(len(short_texts))

with open('short_texts.json', 'w') as f:
    json.dump(short_texts, f)
f.close()

"""Test"""

f = open('short_texts.json')
read_short_texts = json.load(f)
f.close()

print(len(read_short_texts))

s = " "
flag = 0
count = 0
err_op = 0

for i in range(len(read_short_texts)):
    if s in read_short_texts[i]:
        flag = 1
        count += 1

    temp = [x for x in read_short_texts[i] if (x != " ")]
    read_short_texts[i] = temp

    if s in read_short_texts[i]:
        print(i)
        err_op += 1
    else:
        flag = 0 

print(count)
print(err_op)

with open('short_texts.json', 'w') as f:
    json.dump(read_short_texts, f)
f.close()

dict_biterms = {}
dict_terms = {}

for i in range(len(read_short_texts)):
    
    for j in range(len(read_short_texts[i]) - 1):
        
        if read_short_texts[i][j] in dict_terms:
            dict_terms[read_short_texts[i][j]] = dict_terms[read_short_texts[i][j]] + 1
        else:
            dict_terms[read_short_texts[i][j]] = 1
        
        for k in range(j+1, len(read_short_texts[i])):
            str = read_short_texts[i][j] + '_' + read_short_texts[i][k]
            
            if str in dict_biterms:
                dict_biterms[str] = dict_biterms[str] + 1
            else:
                dict_biterms[str] = 1

"""Test"""

with open('dict_terms.json', 'w') as f:
    json.dump(dict_terms, f)
f.close()

with open('dict_biterms.json', 'w') as f:
    json.dump(dict_biterms, f)
f.close()

print(len(dict_biterms))
print(len(dict_terms))

"""Sliding window of 10 words"""

sliding_win = 10
dict_biterms_sw10 = {}

for i in range(len(read_short_texts)):
    for j in range(len(read_short_texts[i]) - 1):
        
        if (j+sliding_win) <= len(read_short_texts[i]):
            kRange = j+sliding_win
        else:
            kRange = len(read_short_texts[i])
        
        for k in range(j+1, kRange):
            str = read_short_texts[i][j] + '_' + read_short_texts[i][k]
            if str in dict_biterms_sw10:
                dict_biterms_sw10[str] = dict_biterms_sw10[str] + 1
            else:
                dict_biterms_sw10[str] = 1

print(len(dict_biterms_sw10))

with open('dict_biterms_sw10.json', 'w') as f:
    json.dump(dict_biterms_sw10, f)
f.close()

dict_biterms = {}
dict_biterms_sw10 = {}

f1 = open('dict_biterms.json')
dict_biterms = json.load(f1)
f1.close()


f3 = open('dict_biterms_sw10.json')
dict_biterms_sw10 = json.load(f3)
f3.close()

num_biterms = len(dict_biterms)

num_biterms_sw10 = len(dict_biterms_sw10)

print("Number of BiTerms in External Dataset: ", num_biterms)
print("Number of BiTerms in External Dataset for sliding window = 10: ", num_biterms_sw10)

dict_terms = {}

f1 = open('dict_terms.json')
dict_terms = json.load(f1)
f1.close()

num_terms = len(dict_terms)
print("Number of Terms in External Dataset: ", num_terms)

total_terms = sum(dict_terms.values())
print("Total terms in External corpus", total_terms)

total_biterms = sum(dict_biterms.values())
print("Total # of Biterms in External corpus", total_biterms)

total_biterms_sw10 = sum(dict_biterms_sw10.values())
print("Total # of Biterms in External corpus with sliding window =10: ", total_biterms_sw10)

pdt = {} #probability_dict_terms

for key, val in dict_terms.items():
    pdt[key] = val / total_terms

print(len(pdt))

pdb = {} #probability_dict_biterms

for key, val in dict_biterms.items():
    pdb[key] = val / total_biterms

print(len(pdb))

pdbsw10 = {} #probability_dict_biterms_sliding_window_10

for key, val in dict_biterms_sw10.items():
    pdbsw10[key] = val / total_biterms_sw10

print(len(pdbsw10))

with open('probability_dict_terms.json', 'w') as f:
    json.dump(pdt, f)
f.close()

with open('probability_dict_biterms.json', 'w') as f:
    json.dump(pdb, f)
f.close()

with open('probability_dict_biterms_sliding_window_10.json', 'w') as f:
    json.dump(pdbsw10, f)
f.close()