""""

              _......_
           .-'      _.-`.
         .'    .-. '-._.-`.
        / /\   |  \        `-..
       / / |   `-.'      .-.   `-.
      /  `-'            (   `.    `.
     |           /\      `-._/      \
     '    .'\   /  `.           _.-'|
    /    /  /   \_.-'        _.':;:/
  .'     \_/             _.-':;_.-'
 /   .-.             _.-' \;.-'
/   (   \       _..-'     |
\    `._/  _..-'          |
 `-.....-'/               |
          |    CUERPOS    \  (o)
     (o)  |                | (\'/)
    (\'/)/                 \;:;
     :;  |                 /)
      ;: `-.._    /__..--'\.' ;:
          :;  `--' :;   :;
"""

class cuerpo:
    def __init__(self, masa, radio):
        self.masa=masa
        self.radio=radio

G=6.67191e-20 #km3 kg-1 s-2

iepetus=cuerpo(1.973e21,734.5)
saturno=cuerpo(5.688e26,58232)
sat=cuerpo(1500, 0.005)

Mi=iepetus.masa #Masa Iepetus [kg]
Ri=iepetus.radio #Radio Iepetus [km]
Ms=saturno.masa #Masa Saturno [kg]
Rs=saturno.radio #Radio Saturno [km]
Dis=3561300 #Distancia Iepetus-Saturno [km]
m=sat.masa #Masa Sat [kg]
Rsoi=Dis*(Mi/Ms)**(2/5) #23317.154948705997 [km]
mu=G*Mi
print(Rsoi)