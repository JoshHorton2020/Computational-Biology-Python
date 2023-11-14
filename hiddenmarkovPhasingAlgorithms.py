import numpy as np
import pandas as pd
import scipy
from plotnine import *


def plot_timeseries_data(data):
    scatter_plot = ggplot(df, aes(x='time', y='level')) + geom_point()
    scatter_plot.save(filename='plot_timeseries_data.png')
    return scatter_plot

approx_duration = 50

'''
Hidden Markov Model maximum likelihood state sequence with 4 states
state 1 - T corresponds to a normal distribution with mean 100 and sd 15
state 2 - A corresponds to a normal dist with mean 150 and sd 25
state 3 - G correcponds to a normal dist with mean 300 and sd 50
state 4 - C corresponds to a normal dist with mean 350 and sd 25
transitions between states are 1/50 and transitions to same state is 49/50
'''

transition_prob = 1/50
unchanging_prob = 49/50

def HMM_MLE(df):
    #inital prob of starting in each state, we have no data so equal prob for each 
    hmmMatrix = [[0, 0, 0, 0]]
    matrixIndex = 0
    states = [scipy.stats.norm(loc=100, scale=15), 
              scipy.stats.norm(loc=150, scale=25), 
              scipy.stats.norm(loc=300, scale=50), 
              scipy.stats.norm(loc=350, scale=25)]
    # x -> row['time'], y -> row['level']
    for index, row in df.iterrows(): 
        currProbs = calculateProbList(hmmMatrix, matrixIndex, states, row)
        hmmMatrix.append(currProbs)
        matrixIndex += 1
    #backtrace
    currY = hmmMatrix[-1].index(max(hmmMatrix[-1]))
    currProb = max(hmmMatrix[-1])
    backtraceMatrix = [currY]
    for col in reversed(hmmMatrix[:-1]):
        tempList = []
        for index, prob in enumerate(col): 
            #you have currY, so calculate the max of each values mult by either trans or unch if index and Y match
            if index is currY: 
                tempList.append(prob + np.log(unchanging_prob))
            else: 
                tempList.append(prob + np.log(transition_prob))
        currProb = max(tempList)
        currY = tempList.index(max(tempList))
        backtraceMatrix = [currY] + backtraceMatrix
    return backtraceMatrix

def calculateProbList(hmmMatrix, matrixIndex, states, row): 
    #gives me probs from previous col
    priorsList = hmmMatrix[matrixIndex]
    #take max mult by prob of data at time given
    maxList = []
    for index, state in enumerate(states): 
        sampleList = []
        for innerdex, prior in enumerate(priorsList): 
            if innerdex is index:
                sampleList.append(prior + np.log(unchanging_prob))
            else: 
                sampleList.append(prior + np.log(transition_prob))
        sampleMax = max(sampleList)
        maxList.append(sampleMax) #-> [tMax, uMax ... ]
    outList = []
    for index, item in enumerate(maxList): 
        #print(states[index].logpdf(row['level'])
        item = item + states[index].logpdf(row['level'])
        outList.append(item)
    return outList
        

def plot_MLE(state_sequence):
    tList = []
    aList = []
    gList = []
    cList = []
    for index, state in enumerate(state_sequence): 
        if state == 0: 
            tList.append([index, state+1])
        elif state == 1: 
            aList.append([index, state+1])
        elif state == 2: 
            gList.append([index, state+1])
        else: 
            cList.append([index, state+1])
    tDf = pd.DataFrame(tList, columns=['index', 'state'])
    aDf = pd.DataFrame(aList, columns=['index', 'state'])
    gDf = pd.DataFrame(gList, columns=['index', 'state'])
    cDf = pd.DataFrame(cList, columns=['index', 'state'])
    plot = ggplot() + geom_point(aes(x='index', y='state'), data=tDf, color='blue') + geom_point(aes(x='index', y='state'), data=aDf, color='red') + geom_point(aes(x='index', y='state'), data=gDf, color='black') + geom_point(aes(x='index', y='state'), data=cDf, color='green')
    plot.save(filename='plot_MLE.png')
    return plot
    
'''
The most likely sequence this data corresponds to given the likely 
event length found from plotting the data
(printing this sequence of A/C/G/Ts)
'''
def MLE_seq(state_sequence, event_length):
    str = ''
    counter = 0
    for index, item in enumerate(state_sequence): 
        if index == 0: 
            counter+=1
            continue
        if item == state_sequence[index-1]:
            counter+=1
        else: 
            numLetters = int(counter/event_length)
            if numLetters == 0: 
                numLetters = 1
            str = str + numToLetter(state_sequence[index-1])*numLetters
            counter = 1
    return str
            
def numToLetter(x): 
    if x == 0: 
        return 'T'
    elif x == 1: 
        return 'A'
    elif x ==2: 
        return 'G'
    else: 
        return 'C'
'''
Forward/backward algorithm giving posterior probabilities for each time point for each level
'''
def HMM_posterior(df):
    #forward
    hmmMatrixForward = [[0, 0, 0, 0]]
    matrixIndex = 0
    states = [scipy.stats.norm(loc=100, scale=15), 
              scipy.stats.norm(loc=150, scale=25), 
              scipy.stats.norm(loc=300, scale=50), 
              scipy.stats.norm(loc=350, scale=25)]
    # x -> row['time'], y -> row['level']
    #going forwward
    for index, row in df.iterrows(): 
        currProbs = forward_backward(hmmMatrixForward, matrixIndex, states, row)
        hmmMatrixForward.append(currProbs)
        matrixIndex += 1
        
    hmmMatrixBackwards = [[1,1,1,1]]
    matrixIndex = 0 
    #going backwards
    for index in range(len(df)-1, -1, -1):
        row = df.iloc[index]
        currProbs = forward_backward(hmmMatrixBackwards, matrixIndex, states, row)
        hmmMatrixBackwards.append(currProbs)
        matrixIndex += 1
    hmmMatrixBackwards = hmmMatrixBackwards[::-1]
    outputList = []
    for outerdex, col in enumerate(hmmMatrixForward): 
        l = []
        for index, probs in enumerate(col): 
            l.append(hmmMatrixForward[outerdex][index] + hmmMatrixBackwards[outerdex][index])
        outputList.append(l)
    return outputList

def forward_backward(hmmMatrix, matrixIndex, states, row): 
    #gives me probs from previous col
    priorsList = hmmMatrix[matrixIndex]
    sumList = []
    for index, state in enumerate(states): 
        sampleList = []
        for innerdex, prior in enumerate(priorsList): 
            if innerdex is index:
                sampleList.append(prior + np.log(unchanging_prob))
            else: 
                sampleList.append(prior + np.log(transition_prob))
        summedVal = scipy.special.logsumexp(sampleList)
        sumList.append(summedVal) #-> [tsum, usum ... ]
    outList = []
    for index, item in enumerate(sumList): 
        #print(states[index].logpdf(row['level'])
        item = item + states[index].logpdf(row['level'])
        outList.append(item)
    return outList

'''
plotting output of posterior probability calculated by our hidden markov model
'''
def plot_posterior(df, posteriors):
    #extract levels from df 
    levels = [row[1]['level'] for row in df.iterrows()]
    #evening the lens
    levels = [0] + levels
    tPost = [subList[0] for subList in posteriors]
    aPost = [subList[1] for subList in posteriors]
    gPost = [subList[2] for subList in posteriors]
    cPost = [subList[3] for subList in posteriors]
    
    data = {'levels': levels, 'T': tPost, 'A': aPost, 'G': gPost, 'C':cPost}
    bigDataFrame = pd.DataFrame(data)
    bigDataFrame = pd.melt(bigDataFrame, id_vars=['levels'], var_name='variable', value_name='value')
    plot = ggplot(bigDataFrame, aes(x='levels', y='value', color='variable')) + facet_wrap('variable') + geom_line()+ xlab('Level') + ylab('Posterior Probability')
    plot.save(filename='plot_posteriors.png')
    return plot

df = pd.read_csv("nanopore.csv")
plot_timeseries_data(df)
state_sequence = HMM_MLE(df)
plot_MLE(state_sequence)
MLE_seq(state_sequence, approx_duration)
posteriors = HMM_posterior(df)
plot_posterior(df, posteriors)

