import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
from playsound import playsound as play
import IPython
from scipy.signal import spectrogram, lfilter, freqz, tf2zpk
from scipy.io import wavfile
#from scipy.signal import spectrogram, lfilter, freqz, tf2zpk{}





s, fs = sf.read('music.wav')
s = s[:250000]
t = np.arange(s.size) / fs
N = s.size
plt.figure(figsize=(6,3))
plt.plot(t, s)

#begin task1
print("delka v s: ",s.size/fs)
print("delka ve vzorcich: ",N)
Max=s[0]
AbsMax=s[0]
Avg=0
for i in range (N):
    Avg+=s[i]
    if AbsMax < abs(s[i]):
        AbsMax= abs(s[i])
    if Max < s[i]:
        Max=s[i]
Avg=Avg/N
print("max: ",Max)
print("Average: ",Avg)
#print("max(|s|): ",AbsMax )


plt.gca().set_xlabel('$t[s]$')
plt.gca().set_title('Zvukový signál')

plt.tight_layout()
plt.show()
#end task1

#begin task2

s-=Avg
s/=AbsMax
#begin sequencing
seq = np.empty((95,1024))
print(len(seq))
arr_cnt=0
pos_cnt=0
cnt=0
i=0
while i<N:
    cnt+=1
    seq[arr_cnt,pos_cnt]=s[i]
    pos_cnt+=1
    if pos_cnt == 1024:
        arr_cnt+=1
        pos_cnt=0
        i=i-512
    i+=1
#end sequencing
#begin plot
sq=90
plt.plot((np.arange(1024) / fs),seq[sq])
plt.show()
seq=seq.T #transpozice
#end polot
#end task2
#begin task3
def fft(s):
    N=len(s)
    if N == 1:
        return s
    M= int(N/2)
    Xeven=[]
    Xodd=[]
    for i in range(M):
        Xeven.append(s[2*i])
        Xodd.append(s[(2*i)+1])
    Feven=fft(Xeven)
    Fodd=fft(Xodd)
    freqbins = np.empty(N, dtype=object)
    for k in range(M):
        #iexp = (-2*np.pi*k/N)*Fodd[k]
        iexp = Fodd[k] * np.exp(1j*(-2*np.pi*k/N))
        freqbins[k]=Feven[k]+iexp
        freqbins[k+M] = Feven[k]-iexp
    return freqbins
odkud = 0     # začátek segmentu v sekundách
odkud_vzorky = int(odkud * fs)         # začátek segmentu ve vzorcích
pokud_vzorky = int(odkud * fs + 1024) # konec segmentu ve vzorcích

s_seg = s[odkud_vzorky:pokud_vzorky]
s_seg_spec = np.fft.fft(s_seg)


G = (1/N * np.abs(s_seg_spec)**2)

_, ax = plt.subplots(3,1)

# np.arange(n) vytváří pole 0..n-1 podobně jako obyč Pythonovský range
ax[0].plot(np.arange(s_seg.size) / fs + odkud, s_seg)
ax[0].set_xlabel('$t[s]$')
ax[0].set_title('Segment signalu $s$')
ax[0].grid(alpha=0.5, linestyle='--')

f = np.arange(G.size) / N * fs
# zobrazujeme prvni pulku spektra
ax[1].plot(range(1024), s_seg_spec)
ax[1].set_xlabel('$f[Hz]$')
ax[1].set_title('Spektralni hustota vykonu [dB]')
ax[1].grid(alpha=0.5, linestyle='--')

s_seg_spec = fft(s_seg)
G = (np.abs(s_seg_spec))

ax[2].plot(range(G.size//2+1), G[:G.size//2+1])
ax[2].set_xlabel('$f[Hz]$')
ax[2].set_title('Spektralni hustota vykonu [dB]')
ax[2].grid(alpha=0.5, linestyle='--')


plt.tight_layout()
plt.show()
#end task3
#begin task4
f, t, sgr = spectrogram(s, fs)
# prevod na PSD
# (ve spektrogramu se obcas objevuji nuly, ktere se nelibi logaritmu, proto +1e-20)
sgr_log = 10 * np.log10(sgr+1e-20)
plt.figure(figsize=(9,3))
plt.pcolormesh(t,f,sgr_log)
plt.gca().set_xlabel('Čas [s]')
plt.gca().set_ylabel('Frekvence [Hz]')
cbar = plt.colorbar()
cbar.set_label('Spektralní hustota výkonu [dB]', rotation=270, labelpad=15)

plt.tight_layout()
plt.show()
#end taks4
#begin task5
#freq "ručně" ...
Gavg =0
cnt=0
cos4= np.zeros((4,2))
for i in range(1024):
    Gavg+= abs(G[i])
Gavg=Gavg/1024
for i in G[:G.size//2+1]:
    if abs(i)>15:
        print(i,(np.where(G==i)[0]*16000/1024)[0])
        cos4[cnt,0]=i
        cos4[cnt,1]=(np.where(G==i)[0]*16000/1024)[0]
        cnt+=1
print(cos4)
"""
fftmax=np.zeros((4,2),dtype=object)
fftmax[0,0]= int(np.where(s_seg_spec==max(s_seg_spec))[0])
fftmax[0,1]= max(s_seg_spec)
s_seg_spec[int(fftmax[0,0])]=0;
s_seg_spec[1024-int(fftmax[0,0])]=0;
fftmax[1,0]= int(np.where(s_seg_spec==max(s_seg_spec))[0])
fftmax[1,1]= max(s_seg_spec)
s_seg_spec[int(fftmax[1,0])]=0;
s_seg_spec[1024-int(fftmax[1,0])]=0;
fftmax[2,0]= int(np.where(s_seg_spec==max(s_seg_spec))[0])
fftmax[2,1]= max(s_seg_spec)
s_seg_spec[int(fftmax[2,0])]=0;
s_seg_spec[int(fftmax[2,0])]=0;
fftmax[3,0]= int(np.where(s_seg_spec==max(s_seg_spec))[0])
fftmax[3,1]= max(s_seg_spec)
s_seg_spec[int(fftmax[3,0])]=0;
s_seg_spec[1024-int(fftmax[3,0])]=0;

print(fftmax)
"""

"""
gde = np.zeros((4))
fftmax = np.sort(s_seg_spec)[0:4]
for i in range(4):
    gde[i]=np.where(s_seg_spec == fftmax[i])
print(fftmax[0])
print(gde)
"""
"""
fftmax = np.zeros((4,2),dtype=complex)#freq,mag
for i in range(1024):
    s_seg_spec_abs=abs(s_seg_spec[i])
    if s_seg_spec_abs > abs(fftmax[0,1]):
        if s_seg_spec_abs >= abs(fftmax[1,1]):
            if s_seg_spec_abs >= abs(fftmax[2,1]): #>1 < 2
                if s_seg_spec_abs >= abs(fftmax[3,1]): #>3
                    np.roll(fftmax,1)
                    fftmax[3]=s_seg_spec[i]
                else: #>2 <3
                    fftmax[0]=fftmax[3]
                    np.roll(fftmax,-1)
                    fftmax[2]=s_seg_spec[i]
            else:
                fftmax[0]=fftmax[1]
                fftmax[1]=s_seg_spec[i]
        else:
            fftmax[0]=s_seg_spec[i]
"""
#end taks5
#begin task6
rat=(16000/1024)
sins = np.zeros((16000*3))
for i in range(len(sins)):
    t=i/16000
    sins[i] = np.sin(2 * np.pi *cos4[0,1] * t ) + np.sin(2 * np.pi *cos4[1,1] * t ) +np.sin(2 * np.pi *cos4[2,1] * t )+np.sin(2 * np.pi *cos4[3,1] * t )
    #cosines and mages and vir done
plt.plot(range(len(sins)),sins)
plt.show()
sf.write('new_file.wav', sins, fs)

s_seg= sins[1024:2048]
s_seg_spec = np.fft.fft(s_seg)
G = (np.abs(s_seg_spec))

plt.plot(range(G.size//2+1), G[:G.size//2+1])
plt.show()
# filtr timew
