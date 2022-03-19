import matplotlib.pyplot as plt
import pandas as pd;

t_2 = [0.181, 0.437, 0.810, 1.219, 1.777, 2.352, 2.795, 3.234, 4.318]
t_4 = [0.171, 0.378, 0.726, 1.026, 1.446, 1.986, 2.625, 3.177, 3.905]
N = [4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000]
#times_N2 = [0.034, 0.072, 0.149, 0.220, 0.587, 1.032, 2.326, 6.438, 12.604, 20.985]
#N_s_N2 = [300, 500, 700, 900, 1500, 2000, 3000, 5000, 7000, 9000]

#data = [[N, time] for N, time in zip(N_s_barnes, times_barnes)]
#data_N2 = [[N, time] for N, time in zip(N_s_N2, times_N2)]
#df1 = pd.DataFrame(data, columns = ['Iterations', 'Error'])
#df2 = pd.DataFrame(data_N2, columns = ['N', 'Time (s)'])
#ax = df1.plot(x='N', y='Time (s)', kind='line', title='Time as a function of N, 100 time steps', \
#     ylabel='Time (s)', legend=True)
#df2.plot(x='N', y='Time (s)', kind='line', title='Time as a function of N, 100 time steps', \
#     ylabel='Time (s)', legend=True)

plt.plot(N, t_2, label='No unrolling')
plt.plot(N, t_4, label='Unrolling')
plt.title("Time vs. equations in the system, 30 iterations")
plt.xlabel('# equations')
plt.ylabel('Time (s)')
plt.legend()
plt.show()
