subroutine heatlut(i,r,g,b)
implicit none
integer :: i
real :: r,g,b
real, dimension(3,256) :: heat = reshape((/0.00000, 0.00000, 0.00000, &
0.01176, 0.00392, 0.00000, &
0.02353, 0.00784, 0.00000, &
0.03529, 0.01176, 0.00000, &
0.04706, 0.01569, 0.00000, &
0.05882, 0.01961, 0.00000, &
0.07059, 0.02353, 0.00000, &
0.08235, 0.02745, 0.00000, &
0.09412, 0.03137, 0.00000, &
0.10588, 0.03529, 0.00000, &
0.11765, 0.03922, 0.00000, &
0.12941, 0.04314, 0.00000, &
0.14118, 0.04706, 0.00000, &
0.15294, 0.05098, 0.00000, &
0.16471, 0.05490, 0.00000, &
0.17647, 0.05882, 0.00000, &
0.18824, 0.06275, 0.00000, &
0.20000, 0.06667, 0.00000, &
0.21176, 0.07059, 0.00000, &
0.22353, 0.07451, 0.00000, &
0.23529, 0.07843, 0.00000, &
0.24706, 0.08235, 0.00000, &
0.25882, 0.08627, 0.00000, &
0.27059, 0.09020, 0.00000, &
0.28235, 0.09412, 0.00000, &
0.29412, 0.09804, 0.00000, &
0.30588, 0.10196, 0.00000, &
0.31765, 0.10588, 0.00000, &
0.32941, 0.10980, 0.00000, &
0.34118, 0.11373, 0.00000, &
0.35294, 0.11765, 0.00000, &
0.36471, 0.12157, 0.00000, &
0.37647, 0.12549, 0.00000, &
0.38824, 0.12941, 0.00000, &
0.40000, 0.13333, 0.00000, &
0.41176, 0.13725, 0.00000, &
0.42353, 0.14118, 0.00000, &
0.43529, 0.14510, 0.00000, &
0.44706, 0.14902, 0.00000, &
0.45882, 0.15294, 0.00000, &
0.47059, 0.15686, 0.00000, &
0.48235, 0.16078, 0.00000, &
0.49412, 0.16471, 0.00000, &
0.50588, 0.16863, 0.00000, &
0.51765, 0.17255, 0.00000, &
0.52941, 0.17647, 0.00000, &
0.54118, 0.18039, 0.00000, &
0.55294, 0.18431, 0.00000, &
0.56471, 0.18824, 0.00000, &
0.57647, 0.19216, 0.00000, &
0.58824, 0.19608, 0.00000, &
0.60000, 0.20000, 0.00000, &
0.61176, 0.20392, 0.00000, &
0.62353, 0.20784, 0.00000, &
0.63529, 0.21176, 0.00000, &
0.64706, 0.21569, 0.00000, &
0.65882, 0.21961, 0.00000, &
0.67059, 0.22353, 0.00000, &
0.68235, 0.22745, 0.00000, &
0.69412, 0.23137, 0.00000, &
0.70588, 0.23529, 0.00000, &
0.71765, 0.23922, 0.00000, &
0.72941, 0.24314, 0.00000, &
0.74118, 0.24706, 0.00000, &
0.75294, 0.25098, 0.00000, &
0.76471, 0.25490, 0.00000, &
0.77647, 0.25882, 0.00000, &
0.78824, 0.26275, 0.00000, &
0.80000, 0.26667, 0.00000, &
0.81176, 0.27059, 0.00000, &
0.82353, 0.27451, 0.00000, &
0.83529, 0.27843, 0.00000, &
0.84706, 0.28235, 0.00000, &
0.85882, 0.28627, 0.00000, &
0.87059, 0.29020, 0.00000, &
0.88235, 0.29412, 0.00000, &
0.89412, 0.29804, 0.00000, &
0.90588, 0.30196, 0.00000, &
0.91765, 0.30588, 0.00000, &
0.92941, 0.30980, 0.00000, &
0.94118, 0.31373, 0.00000, &
0.95294, 0.31765, 0.00000, &
0.96471, 0.32157, 0.00000, &
0.97647, 0.32549, 0.00000, &
0.98824, 0.32941, 0.00000, &
1.00000, 0.33333, 0.00000, &
1.00000, 0.33725, 0.00000, &
1.00000, 0.34118, 0.00000, &
1.00000, 0.34510, 0.00000, &
1.00000, 0.34902, 0.00000, &
1.00000, 0.35294, 0.00000, &
1.00000, 0.35686, 0.00000, &
1.00000, 0.36078, 0.00000, &
1.00000, 0.36471, 0.00000, &
1.00000, 0.36863, 0.00000, &
1.00000, 0.37255, 0.00000, &
1.00000, 0.37647, 0.00000, &
1.00000, 0.38039, 0.00000, &
1.00000, 0.38431, 0.00000, &
1.00000, 0.38824, 0.00000, &
1.00000, 0.39216, 0.00000, &
1.00000, 0.39608, 0.00000, &
1.00000, 0.40000, 0.00000, &
1.00000, 0.40392, 0.00000, &
1.00000, 0.40784, 0.00000, &
1.00000, 0.41176, 0.00000, &
1.00000, 0.41569, 0.00000, &
1.00000, 0.41961, 0.00000, &
1.00000, 0.42353, 0.00000, &
1.00000, 0.42745, 0.00000, &
1.00000, 0.43137, 0.00000, &
1.00000, 0.43529, 0.00000, &
1.00000, 0.43922, 0.00000, &
1.00000, 0.44314, 0.00000, &
1.00000, 0.44706, 0.00000, &
1.00000, 0.45098, 0.00000, &
1.00000, 0.45490, 0.00000, &
1.00000, 0.45882, 0.00000, &
1.00000, 0.46275, 0.00000, &
1.00000, 0.46667, 0.00000, &
1.00000, 0.47059, 0.00000, &
1.00000, 0.47451, 0.00000, &
1.00000, 0.47843, 0.00000, &
1.00000, 0.48235, 0.00000, &
1.00000, 0.48627, 0.00000, &
1.00000, 0.49020, 0.00000, &
1.00000, 0.49412, 0.00000, &
1.00000, 0.49804, 0.00000, &
1.00000, 0.50196, 0.00000, &
1.00000, 0.50588, 0.00000, &
1.00000, 0.50980, 0.00000, &
1.00000, 0.51373, 0.00000, &
1.00000, 0.51765, 0.00000, &
1.00000, 0.52157, 0.00000, &
1.00000, 0.52549, 0.00000, &
1.00000, 0.52941, 0.00000, &
1.00000, 0.53333, 0.00000, &
1.00000, 0.53725, 0.00000, &
1.00000, 0.54118, 0.00000, &
1.00000, 0.54510, 0.00000, &
1.00000, 0.54902, 0.00000, &
1.00000, 0.55294, 0.00000, &
1.00000, 0.55686, 0.00000, &
1.00000, 0.56078, 0.00000, &
1.00000, 0.56471, 0.00000, &
1.00000, 0.56863, 0.00000, &
1.00000, 0.57255, 0.00000, &
1.00000, 0.57647, 0.00000, &
1.00000, 0.58039, 0.00000, &
1.00000, 0.58431, 0.00000, &
1.00000, 0.58824, 0.00000, &
1.00000, 0.59216, 0.00000, &
1.00000, 0.59608, 0.00000, &
1.00000, 0.60000, 0.00000, &
1.00000, 0.60392, 0.00000, &
1.00000, 0.60784, 0.00000, &
1.00000, 0.61176, 0.00000, &
1.00000, 0.61569, 0.00000, &
1.00000, 0.61961, 0.00000, &
1.00000, 0.62353, 0.00000, &
1.00000, 0.62745, 0.00000, &
1.00000, 0.63137, 0.00000, &
1.00000, 0.63529, 0.00000, &
1.00000, 0.63922, 0.00000, &
1.00000, 0.64314, 0.00000, &
1.00000, 0.64706, 0.00000, &
1.00000, 0.65098, 0.01176, &
1.00000, 0.65490, 0.02353, &
1.00000, 0.65882, 0.03529, &
1.00000, 0.66275, 0.04706, &
1.00000, 0.66667, 0.05882, &
1.00000, 0.67059, 0.07059, &
1.00000, 0.67451, 0.08235, &
1.00000, 0.67843, 0.09412, &
1.00000, 0.68235, 0.10588, &
1.00000, 0.68627, 0.11765, &
1.00000, 0.69020, 0.12941, &
1.00000, 0.69412, 0.14118, &
1.00000, 0.69804, 0.15294, &
1.00000, 0.70196, 0.16471, &
1.00000, 0.70588, 0.17647, &
1.00000, 0.70980, 0.18824, &
1.00000, 0.71373, 0.20000, &
1.00000, 0.71765, 0.21176, &
1.00000, 0.72157, 0.22353, &
1.00000, 0.72549, 0.23529, &
1.00000, 0.72941, 0.24706, &
1.00000, 0.73333, 0.25882, &
1.00000, 0.73725, 0.27059, &
1.00000, 0.74118, 0.28235, &
1.00000, 0.74510, 0.29412, &
1.00000, 0.74902, 0.30588, &
1.00000, 0.75294, 0.31765, &
1.00000, 0.75686, 0.32941, &
1.00000, 0.76078, 0.34118, &
1.00000, 0.76471, 0.35294, &
1.00000, 0.76863, 0.36471, &
1.00000, 0.77255, 0.37647, &
1.00000, 0.77647, 0.38824, &
1.00000, 0.78039, 0.40000, &
1.00000, 0.78431, 0.41176, &
1.00000, 0.78824, 0.42353, &
1.00000, 0.79216, 0.43529, &
1.00000, 0.79608, 0.44706, &
1.00000, 0.80000, 0.45882, &
1.00000, 0.80392, 0.47059, &
1.00000, 0.80784, 0.48235, &
1.00000, 0.81176, 0.49412, &
1.00000, 0.81569, 0.50588, &
1.00000, 0.81961, 0.51765, &
1.00000, 0.82353, 0.52941, &
1.00000, 0.82745, 0.54118, &
1.00000, 0.83137, 0.55294, &
1.00000, 0.83529, 0.56471, &
1.00000, 0.83922, 0.57647, &
1.00000, 0.84314, 0.58824, &
1.00000, 0.84706, 0.60000, &
1.00000, 0.85098, 0.61176, &
1.00000, 0.85490, 0.62353, &
1.00000, 0.85882, 0.63529, &
1.00000, 0.86275, 0.64706, &
1.00000, 0.86667, 0.65882, &
1.00000, 0.87059, 0.67059, &
1.00000, 0.87451, 0.68235, &
1.00000, 0.87843, 0.69412, &
1.00000, 0.88235, 0.70588, &
1.00000, 0.88627, 0.71765, &
1.00000, 0.89020, 0.72941, &
1.00000, 0.89412, 0.74118, &
1.00000, 0.89804, 0.75294, &
1.00000, 0.90196, 0.76471, &
1.00000, 0.90588, 0.77647, &
1.00000, 0.90980, 0.78824, &
1.00000, 0.91373, 0.80000, &
1.00000, 0.91765, 0.81176, &
1.00000, 0.92157, 0.82353, &
1.00000, 0.92549, 0.83529, &
1.00000, 0.92941, 0.84706, &
1.00000, 0.93333, 0.85882, &
1.00000, 0.93725, 0.87059, &
1.00000, 0.94118, 0.88235, &
1.00000, 0.94510, 0.89412, &
1.00000, 0.94902, 0.90588, &
1.00000, 0.95294, 0.91765, &
1.00000, 0.95686, 0.92941, &
1.00000, 0.96078, 0.94118, &
1.00000, 0.96471, 0.95294, &
1.00000, 0.96863, 0.96471, &
1.00000, 0.97255, 0.97647, &
1.00000, 0.97647, 0.98824, &
1.00000, 0.98039, 1.00000, &
1.00000, 0.98431, 1.00000, &
1.00000, 0.98824, 1.00000, &
1.00000, 0.99216, 1.00000, &
1.00000, 0.99608, 1.00000, &
1.00000, 1.00000, 1.00000 /),(/3,256/))

if(i.le.0)then
   r=0.0
   g=0.0
   b=0.0
elseif(i.gt.256)then
   r=1.0
   g=1.0
   b=1.0
else
   r=heat(1,i)
   g=heat(2,i)
   b=heat(3,i)
endif

end subroutine heatlut
