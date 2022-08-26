# Numerical-analysis-example-GUI.
GUI visualisation of plotting funtion derivative, calculated using Runge-Kutta and Adams methods.

A sandbox project for getting familiar with **tkinter**.

A numerical solution of the differential equation $\dfrac{dz_1(t)}{dt}=2 \cdot z_1(t) - t^2$, boarder condtition $z_1(0)=0.25$ is calculated by Runge-Kutta and Adams methods. A user can change time grid step $h$ and $b$, end of interval $(a=0,b)$ of the $t$ axis.The resulting derivative is calculated. Then the area between the resulting function and the function $z_2(t)=0.25+5 \cdot \ln(t+1)^2 \cdot \dfrac{\sin(0.5t)}{0.1}$ is computed using trapezoidal rule. Also the points of intersection of two funcions are found.

The GUI was made using **tkinter**.
