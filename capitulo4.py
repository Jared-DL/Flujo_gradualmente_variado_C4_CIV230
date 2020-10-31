def inicio():
    print(""" BIENVENIDO A CÁLCULO DE FLUJO GRADUALMENTE VARIADO
PARA SALIR EN CUALQUIER MOMENTO ESCRIBE: -EXIT()-
---*---*---*---Jared DL---*---*---*---""")
    seleccion_de_seccion()

def seleccion_de_seccion() :
    global seccion_transversal
    seccion_transversal = int(input("""
Tipo de sección transversal:
1. Rectangular
2. Trapecial
Seleciona el tipo de sección transversal: """))
    if seccion_transversal == 1:
        print('SELECCIONASTE SECCIÓN RETANGULAR')
        datos_principales()
    elif seccion_transversal == 2:
        print('SELECCIONASTE SECCIÓN TRAPECIAL')
        datos_principales()
    else :
        print('SELECCIONA UNA OPCIÓN VÁLIDA!!!')
        seleccion_de_seccion()

    

def datos_principales() :
    global Q_caudal
    global g_gravedad
    global n_coeficiente_de_manning
    global So_pendiente_solera
    global Y_dato
    print("Nesecitamos los siguientes datos:")
    Q_caudal=float(input("Caudal Q [m3/s]= "))
    g_gravedad=float(input("Gravedad [m/s2]= "))
    n_coeficiente_de_manning=float(input("Coeficiente de manning n = "))
    So_pendiente_solera=float(input("Pendiente de la solera So [m/m]= "))
    Y_dato = float(input("Y tirante dato = "))
    if seccion_transversal==1 :
        rectangular_ingresar_datos()
    elif seccion_transversal==2 :
        trapecial_ingresar_datos()


def rectangular_ingresar_datos() :
    print("Datos para la seccion rectangular")
    global b_solera
    b_solera=float(input("Solera b[m] = "))
    calculos_seccion_rectangular()


def trapecial_ingresar_datos() :
    print("Datos para la seccion trapecial")
    global b_solera
    global z_lateral_trapecio
    b_solera=float(input("Solera b[m] = "))
    z_lateral_trapecio=float(input("Ancho por metro de altura de las paredes laterales Z[m] = "))
    calculos_seccion_trapecial()

def calculos_seccion_rectangular() :
    global Yc_tirante_critico
    global Ac_area_critica
    global Pc_perimetro_critico
    global Sc_pendiente_critica
    global Yn_tirante_normal
    Yc_tirante_critico=(Q_caudal**2/(g_gravedad*b_solera**2))**(1/3)
    print(f'Tirante Crítico Yc = {Yc_tirante_critico} [m]')
    Ac_area_critica=calcular_area_rectangulo(b_solera,Yc_tirante_critico)
    print(f'Área crítica Ac = {Ac_area_critica} [m2]')
    Pc_perimetro_critico =calcular_perimetro_mojado_rectangulo(b_solera,Yc_tirante_critico)
    print(f'Perímetro crítico Pc = {Pc_perimetro_critico} [m2]')
    Sc_pendiente_critica = calcular_pendiente(Q_caudal,n_coeficiente_de_manning,Pc_perimetro_critico,Ac_area_critica)
    print(f'Pendiente Crítica Sc= {Sc_pendiente_critica} [m/m]')
    Yn_tirante_normal=newton_rahpson_tirate_normal_para_seccion_rectangular()
    print(f'Tirante normal= {Yn_tirante_normal} [m/m]')    
    determinar_pendiente_del_canal()
    determinar_zona_del_perfil(pendiente_del_canal)
    determinar_tipo_de_perfil()
    eleccion_del_metodo()


def calculos_seccion_trapecial() :
    global Yc_tirante_critico
    global Ac_area_critica
    global Pc_perimetro_critico
    global Sc_pendiente_critica
    global Yn_tirante_normal
    Yc_tirante_critico=newton_rahpson_tirate_critico_para_seccion_trapecial()
    print(f'Tirante Crítico Yc = {Yc_tirante_critico} [m]')
    Ac_area_critica=calcular_area_trapecio(b_solera,Yc_tirante_critico,z_lateral_trapecio)
    print(f'Área crítica Ac = {Ac_area_critica} [m2]')
    Pc_perimetro_critico =calcular_perimetro_mojado_trapecio(b_solera,Yc_tirante_critico,z_lateral_trapecio)
    print(f'Perímetro crítico Pc = {Pc_perimetro_critico} [m2]')
    Sc_pendiente_critica = calcular_pendiente(Q_caudal,n_coeficiente_de_manning,Pc_perimetro_critico,Ac_area_critica)
    print(f'Pendiente Crítica Sc= {Sc_pendiente_critica} [m/m]')
    Yn_tirante_normal=newton_rahpson_tirate_normal_para_seccion_trapecial()
    print(f'Tirante normal= {Yn_tirante_normal} [m/m]')    
    determinar_pendiente_del_canal()
    determinar_zona_del_perfil(pendiente_del_canal)
    determinar_tipo_de_perfil()
    eleccion_del_metodo()

def newton_rahpson_tirate_normal_para_seccion_rectangular() :
    yn=1.1
    yn_1=0.1
    contador=0
    while contador<20 :
        funcion_y=((So_pendiente_solera**(1/2)/n_coeficiente_de_manning)*((b_solera*yn)**(5/3))/(2*yn+b_solera)**(2/3))-Q_caudal
        funcion_y_derivada=(So_pendiente_solera**(1/2)/n_coeficiente_de_manning)*(b_solera*(b_solera*yn)**(2/3)*(6*yn+5*b_solera))/(3*(2*yn+b_solera)**(5/3))
        yn_1=yn-funcion_y/funcion_y_derivada
        yn=yn_1
        contador += 1
    return yn

def newton_rahpson_tirate_critico_para_seccion_trapecial() :
    yc=1
    contador=0
    while contador<20 :
        funcion_y=  ((((b_solera+z_lateral_trapecio*yc)*yc)**(3))/(b_solera+2*z_lateral_trapecio*yc))-((Q_caudal**(2))/(g_gravedad))
        funcion_y_derivada= ((b_solera**(6))/(32*(2*yc*z_lateral_trapecio+b_solera)**(2)*z_lateral_trapecio**(2)))+((5*yc**(4)*z_lateral_trapecio**(2))/(2))+5*b_solera*yc**(3)*z_lateral_trapecio+((21*b_solera**(2)*yc**(2))/(8))+((b_solera**(3)*yc)/(8*z_lateral_trapecio))-((b_solera**(4))/(32*z_lateral_trapecio**(2)))
        yc_1=yc-funcion_y/funcion_y_derivada
        yc=yc_1
        contador += 1
    return yc

def newton_rahpson_tirate_normal_para_seccion_trapecial() :
    yn=1.1
    contador=0
    while contador<20 :
        funcion_y=((So_pendiente_solera**(((1)/(2))))/(n_coeficiente_de_manning))*((((b_solera+z_lateral_trapecio*yn)*yn)**(((5)/(3))))/((2*yn*(z_lateral_trapecio**(2)+1)**(((1)/(2)))+b_solera)**(((2)/(3)))))-Q_caudal
        funcion_y_derivada=(((So_pendiente_solera)**(1/2)*(16*yn**(2)*z_lateral_trapecio*(z_lateral_trapecio**(2)+1)**(1/2)+2*b_solera*yn*(3*(z_lateral_trapecio**(2)+1)**(1/2)+5*z_lateral_trapecio)+5*b_solera**(2))*(yn*(yn*z_lateral_trapecio+b_solera))**(((2)/(3))))/(3*n_coeficiente_de_manning*(2*yn*(z_lateral_trapecio**(2)+1)**(1/2)+b_solera)**(((5)/(3)))))
        yn_1=yn-funcion_y/funcion_y_derivada
        yn=yn_1
        contador += 1
    return yn   


def determinar_pendiente_del_canal() :
    global pendiente_del_canal
    diferencia_de_soleras=abs(Sc_pendiente_critica-So_pendiente_solera)
    print(f'diferencia de soleras: {diferencia_de_soleras}')
    if ((diferencia_de_soleras)<0.001) :
        print(f'{So_pendiente_solera}={Sc_pendiente_critica}  ==> So=Sc Perfil Tipo C critico')
        pendiente_del_canal='C'
    elif So_pendiente_solera<Sc_pendiente_critica :
        print(f'{So_pendiente_solera}<{Sc_pendiente_critica} ==> So<Sc Perfil Tipo M Mild')
        pendiente_del_canal='M'
    elif So_pendiente_solera>Sc_pendiente_critica :
        print(f'{So_pendiente_solera}>{Sc_pendiente_critica} ==> So>Sc Perfil Tipo S Strong')
        pendiente_del_canal='S'
    elif So_pendiente_solera==0 :
        print(f'{So_pendiente_solera}=So Perfil Tipo H Horizontal')
        pendiente_del_canal='H'
    elif So_pendiente_solera<0 :
        print(f'{So_pendiente_solera}=So ==> So<o Perfil Tipo A Adverse')
        pendiente_del_canal='A'
    else :
        print('Algo anda mal, el canal no encaja con ningún tipo conocido')

def determinar_zona_del_perfil(pendiente_del_canal) :
    global zona_del_perfil
    if pendiente_del_canal=='C' :
        if Y_dato>Yc_tirante_critico :
            zona_del_perfil='1'
            print(f'La zona del perfil es: {zona_del_perfil}')
        elif abs(Yc_tirante_critico-Y_dato)<0.001:
            zona_del_perfil='2'
            print(f'La zona del perfil es: {zona_del_perfil}')
        elif Y_dato<Yc_tirante_critico:
            zona_del_perfil='3'
            print(f'La zona del perfil es: {zona_del_perfil}')
        else :
            print('ocurrio un error')          
    if ((Y_dato>Yn_tirante_normal) and (Y_dato>Yc_tirante_critico)) :
        zona_del_perfil='1'
        print(f'La zona del perfil es: {zona_del_perfil}')
    elif((Yn_tirante_normal>Y_dato) and (Y_dato>Yc_tirante_critico)) :
        zona_del_perfil='2'
        print(f'La zona del perfil es: {zona_del_perfil}')
    elif ((Yn_tirante_normal>Yc_tirante_critico)and(Yc_tirante_critico>Y_dato)) :
        zona_del_perfil='3'
        print(f'La zona del perfil es: {zona_del_perfil}')
    else :
        print('NO se pudo determinar la zona de perfil ocurrio un error')

def determinar_tipo_de_perfil() :
    tipo_de_perfil = pendiente_del_canal+zona_del_perfil
    global tipo_de_curva
    global tipo_de_flujo
    global sentido_de_calculo
    if tipo_de_perfil == 'M1' or tipo_de_perfil=='M3' or tipo_de_perfil =='C1' or tipo_de_perfil =='C3' or tipo_de_perfil =='S1' or tipo_de_perfil =='S3' or tipo_de_perfil =='H3' or tipo_de_perfil =='A3':
        tipo_de_curva='Remanso'
    elif tipo_de_perfil == 'M2' or tipo_de_perfil =='S2' or tipo_de_perfil =='H2' or tipo_de_perfil =='A2' :
        tipo_de_curva = 'Caída'
    elif tipo_de_perfil == 'C2':
        tipo_de_curva = 'paralelo al fondo'
    elif tipo_de_perfil =='H1' or tipo_de_perfil =='A1':
        tipo_de_curva= 'No existe'
    else:
        tipo_de_curva='No se pudo determinar el tipo de curva Ocurrio un error'
    if tipo_de_perfil == 'M1' or tipo_de_perfil =='C1' or tipo_de_perfil =='S1' or tipo_de_perfil =='H2' or tipo_de_perfil =='A2' or tipo_de_perfil =='M2' :
        tipo_de_flujo = 'Subcrítico'
    elif tipo_de_perfil == 'M3' or tipo_de_perfil =='C3' or tipo_de_perfil =='S2' or tipo_de_perfil =='S3' or tipo_de_perfil =='H3' or tipo_de_perfil =='A3' :
        tipo_de_flujo = 'Supercrítico'
    elif tipo_de_perfil == 'C2':
        tipo_de_flujo = 'Crítico'
    elif tipo_de_perfil == 'A1' or tipo_de_perfil =='H1' :
        tipo_de_flujo = 'no existe' 
    else :
        tipo_de_flujo = 'No sepudo determinar el tipo de flujo ocurrio un error'
    if tipo_de_perfil == 'M1' or tipo_de_perfil =='M2' or tipo_de_perfil =='S1' or tipo_de_perfil =='C1' or tipo_de_perfil =='H2' or tipo_de_perfil =='A2' :
        sentido_de_calculo='De aguas abajo HACIA aguas arriba'
    elif tipo_de_perfil == 'M3' or tipo_de_perfil =='S2' or tipo_de_perfil =="S3" or tipo_de_perfil =='C3' or tipo_de_perfil =='H3' or tipo_de_perfil =='A3' :
        sentido_de_calculo='De aguas arriba HACIA aguas abajo'
    elif tipo_de_perfil == 'C2' :
        sentido_de_calculo ='El tirante Y es constaante'
    elif tipo_de_perfil == 'A1' or tipo_de_perfil =='H1' :
        sentido_de_calculo = 'No se pudo determinar el sentido de calculo ERROR'
    else :
        sentido_de_calculo= "ocurrio un error"
    print(f'El tipo de perfil es: {tipo_de_perfil}')
    print(f'El tipo de curva es: {tipo_de_curva}')
    print(f'El tipo de flujo es: {tipo_de_flujo}')
    print(f'Sentido de calculo es: {sentido_de_calculo}')

def eleccion_del_metodo() :
    global metodo_elegido
    opcion_de_metodo_elegido=int(input(""" Elige qué método de resolución a seguir:
0. Método directo por tramos(Hallar Delta X dado Yo y Yf)
1. Método estándar por pasos(Hallar Yf dado Yo y DeltaX)
...: """))
    if opcion_de_metodo_elegido == 0 :
        metodo_elegido='Directo por tramos'
        if seccion_transversal == 1 :
            calcular_delta_x_sección_rectangular()
        elif seccion_transversal == 2 :
            calcular_delta_x_sección_trapecial()
    elif opcion_de_metodo_elegido == 1 :
        metodo_elegido='Directo por tramos'
        print('RECUERDA QUE SI SENTIDO DE CALCULO ES DE AGUAS ABAJO HACIA AGUAS ARRIBA DELTA X SE DEBE ANOTAR CON SIGNO NEGATIVO')
        if seccion_transversal == 1 :
            calcular_tirante_sección_rectangular()
        elif seccion_transversal == 2 :
            calcular_tirante_sección_trapecial()
    else :
        print('elige una opcionn válida!!!')
        eleccion_del_metodo()

# MÉTODO DIRECTO POR TRAMOS:

def calcular_delta_x_sección_rectangular() :
    print('Hallaremos la distancia X entre dos puntos del canal:')
    print('recuerda que si la dirección de cálculo es: -De aguas abajo HACIA aguas arriba- Y1 debe estar aguas abajo y Y2 debe estar aguas arriba')
    y1=float(input('Valor del tirante 1 Y1[m]='))
    y2=float(input('Valor del tirante 2 Y2[m]='))
    A1=calcular_area_rectangulo(b_solera,y1)
    P1=calcular_perimetro_mojado_rectangulo(b_solera,y1)
    E1=calcular_energia(g_gravedad, y1, Q_caudal, A1)
    S1=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,P1,A1)
    A2=calcular_area_rectangulo(b_solera,y2)
    P2=calcular_perimetro_mojado_rectangulo(b_solera,y2)
    E2=calcular_energia(g_gravedad, y2, Q_caudal, A2)
    S2=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,P2,A2)
    Sf=(S1+S2)/2
    delta_x=(E2-E1)/(So_pendiente_solera-Sf)
    print(f'A1[m2]= {A1}')
    print(f'P1[m]= {P1}')
    print(f'E1[m]= {E1}')
    print(f'S1[m/m]= {S1}')
    print(f'A2[m2]= {A2}')
    print(f'P2[m]= {P2}')
    print(f'E2[m]= {E2}')
    print(f'S2[m/m]= {S2}')
    print(f'Sf[m/m]= {Sf}')
    print(f'Delta x[m]= {delta_x}')
    preguntar_si_se_desea_hacer_mas_calculos_de_deltas_x()

def calcular_delta_x_sección_trapecial() :
    print('Hallaremos la distancia X entre dos puntos del canal:')
    print('recuerda que si la dirección de cálculo es: -De aguas abajo HACIA aguas arriba- Y1 debe estar aguas abajo y Y2 debe estar aguas arriba')
    y1=float(input('Valor del tirante 1 Y1[m]='))
    y2=float(input('Valor del tirante 2 Y2[m]='))
    A1=calcular_area_trapecio(b_solera,y1,z_lateral_trapecio)
    P1=calcular_perimetro_mojado_trapecio(b_solera,y1,z_lateral_trapecio)
    E1=calcular_energia(g_gravedad,y1,Q_caudal,A1)
    S1=calcular_pendiente(Q_caudal, n_coeficiente_de_manning,P1,A1)
    A2=calcular_area_trapecio(b_solera,y2,z_lateral_trapecio)
    P2=calcular_perimetro_mojado_trapecio(b_solera,y2,z_lateral_trapecio)
    E2=calcular_energia(g_gravedad,y2,Q_caudal,A2)
    S2=calcular_pendiente(Q_caudal, n_coeficiente_de_manning,P2,A2)
    Sf=(S1+S2)/2
    delta_x=(E2-E1)/(So_pendiente_solera-Sf)
    print(f'A1[m2]= {A1}')
    print(f'P1[m]= {P1}')
    print(f'E1[m]= {E1}')
    print(f'S1[m/m]= {S1}')
    print(f'A2[m2]= {A2}')
    print(f'P2[m]= {P2}')
    print(f'E2[m]= {E2}')
    print(f'S2[m/m]= {S2}')
    print(f'Sf[m/m]= {Sf}')
    print(f'Delta x[m]= {delta_x}')
    preguntar_si_se_desea_hacer_mas_calculos_de_deltas_x()


def preguntar_si_se_desea_hacer_mas_calculos_de_deltas_x() :
    volver_a_calcular_delta_x = int(input("""Deseas calcular mas deltas?
0.NO (si eliges esta opción el programa finalizará)
1.SI
Elige una opción: """))
    if volver_a_calcular_delta_x == 1 :
        if seccion_transversal == 1 :
            calcular_delta_x_sección_rectangular()
        elif seccion_transversal == 2 :
            calcular_delta_x_sección_trapecial()
    elif volver_a_calcular_delta_x == 0 :
        print('Gracias por usar mi Programa, atte: Jared DL')
    else :
        print('Seleciona una opcion válida!!!')
        preguntar_si_se_desea_hacer_mas_calculos_de_deltas_x()

# ESTÁNDAR POR PASOS:

def calcular_tirante_sección_rectangular() :
    Y1=float(input('Valor del tirante 1 Y1[m]='))
    Z1=float(input('Valor de altura 1 Z1[m]='))
    delta_x=float(input('Valor de Delta X [m]='))
    Z2=float(input('Valor de altura 2 Z2[m]='))
    Y2=calcular_H2_para_seccion_rectangular(Y1,Z1,Z2,delta_x)
    print(f'el tirante´para un DeltaX={delta_x} es Y2={Y2} y el delta H2[m]={delta_H2_final}')
    preguntar_si_se_desea_hacer_mas_calculos_de_tirantes()
    

def calcular_H2_para_seccion_rectangular(tirante_inicial,altura_del_canal_inicial,altura_del_canal_final,delta_x):
    A1=calcular_area_rectangulo(b_solera,tirante_inicial)
    P1=calcular_perimetro_mojado_rectangulo(b_solera,tirante_inicial)
    H1=calcular_altura_de_energía(altura_del_canal_inicial,tirante_inicial,Q_caudal,A1)
    S1=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,P1,A1)
    tirante_final=0.1
    tirante_final=hallar_tirante_final_con_error_rectangular(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.1)
    tirante_final-=0.1
    tirante_final=hallar_tirante_final_con_error_rectangular(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.01)
    tirante_final-=0.01
    tirante_final=hallar_tirante_final_con_error_rectangular(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.001)
    tirante_final-=0.001
    tirante_final=hallar_tirante_final_con_error_rectangular(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.0001)
    tirante_final-=0.0001
    tirante_final=hallar_tirante_final_con_error_rectangular(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.00001)
    A2=calcular_area_rectangulo(b_solera,tirante_final)
    P2=calcular_perimetro_mojado_rectangulo(b_solera,tirante_final)
    H2_1=calcular_altura_de_energía(altura_del_canal_inicial,tirante_final,Q_caudal,A2)
    S2=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,P2,A2)
    Sf=(S1+S2)/2
    Hf=Sf*delta_x
    H2_2=H1-Hf
    print(f'A1[m2]= {A1}')
    print(f'P1[m]= {P1}')
    print(f'H1[m]= {H1}')
    print(f'S1[m/m]= {S1}')
    print(f'A2[m2]= {A2}')
    print(f'P2[m]= {P2}')
    print(f'H2_1[m]= {H2_1}')
    print(f'S2[m/m]= {S2}')
    print(f'Sf[m/m]= {Sf}')
    print(f'Hf[m/m]= {Hf}')
    print(f'H2_2[m]= {H2_2}')
    print(f'Delta x[m]= {delta_x}')
    return tirante_final

def hallar_tirante_final_con_error_rectangular(tirante,pendiente_1,altura_de_energia_1,altura_del_canal_final,delta_x,error) :
    global delta_H2_final
    delta_H2=1
    tirante_final_con_error=tirante
    while delta_H2>error :
        area=calcular_area_rectangulo(b_solera,tirante_final_con_error)
        perimetro=calcular_perimetro_mojado_rectangulo(b_solera,tirante_final_con_error)
        H2_1=calcular_altura_de_energía(altura_del_canal_final,tirante_final_con_error,Q_caudal,area)
        S2=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,perimetro,area)
        Sf=(pendiente_1+S2)/2
        Hf=Sf*delta_x
        H2_2=altura_de_energia_1-Hf
        delta_H2=abs(H2_1-H2_2)
        delta_H2_final=delta_H2
        tirante_final_con_error += error
    return tirante_final_con_error


def calcular_tirante_sección_trapecial() :
    Y1=float(input('Valor del tirante 1 Y1[m]='))
    Z1=float(input('Valor de altura 1 Z1[m]='))
    delta_x=float(input('Valor de Delta X [m]='))
    Z2=float(input('Valor de altura 2 Z2[m]='))
    Y2=calcular_H2_para_seccion_trapecial(Y1,Z1,Z2,delta_x)
    print(f'El tirante para un DeltaX={delta_x} es Y2={Y2} y el delta H2[m]={delta_H2_final}')
    
def calcular_H2_para_seccion_trapecial(tirante_inicial,altura_del_canal_inicial,altura_del_canal_final,delta_x):
    A1=calcular_area_trapecio(b_solera,tirante_inicial,z_lateral_trapecio)
    P1=calcular_perimetro_mojado_trapecio(b_solera,tirante_inicial,z_lateral_trapecio)
    H1=calcular_altura_de_energía(altura_del_canal_inicial,tirante_inicial,Q_caudal,A1)
    S1=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,P1,A1)
    tirante_final=0.1
    tirante_final=hallar_tirante_final_con_error_trapecial(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.1)
    tirante_final-=0.1
    tirante_final=hallar_tirante_final_con_error_trapecial(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.01)
    tirante_final-=0.01
    tirante_final=hallar_tirante_final_con_error_trapecial(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.001)
    tirante_final-=0.001
    tirante_final=hallar_tirante_final_con_error_trapecial(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.0001)
    tirante_final-=0.0001
    tirante_final=hallar_tirante_final_con_error_trapecial(tirante_final,S1,H1,altura_del_canal_final,delta_x,0.00001)
    A2=calcular_area_trapecio(b_solera,tirante_final,z_lateral_trapecio)
    P2=calcular_perimetro_mojado_trapecio(b_solera,tirante_final,z_lateral_trapecio)
    H2_1=calcular_altura_de_energía(altura_del_canal_inicial,tirante_final,Q_caudal,A2)
    S2=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,P2,A2)
    Sf=(S1+S2)/2
    Hf=Sf*delta_x
    H2_2=H1-Hf
    print(f'A1[m2]= {A1}')
    print(f'P1[m]= {P1}')
    print(f'H1[m]= {H1}')
    print(f'S1[m/m]= {S1}')
    print(f'A2[m2]= {A2}')
    print(f'P2[m]= {P2}')
    print(f'H2_1[m]= {H2_1}')
    print(f'S2[m/m]= {S2}')
    print(f'Sf[m/m]= {Sf}')
    print(f'Hf[m/m]= {Hf}')
    print(f'H2_2[m]= {H2_2}')
    print(f'Delta x[m]= {delta_x}')
    return tirante_final

def hallar_tirante_final_con_error_trapecial(tirante,pendiente_1,altura_de_energia_1,altura_del_canal_final,delta_x,error) :
    global delta_H2_final
    delta_H2=1
    tirante_final_con_error=tirante
    while delta_H2>error :
        area=calcular_area_trapecio(b_solera,tirante_final_con_error,z_lateral_trapecio)
        perimetro=calcular_perimetro_mojado_trapecio(b_solera,tirante_final_con_error,z_lateral_trapecio)
        H2_1=calcular_altura_de_energía(altura_del_canal_final,tirante_final_con_error,Q_caudal,area)
        S2=calcular_pendiente(Q_caudal,n_coeficiente_de_manning,perimetro,area)
        Sf=(pendiente_1+S2)/2
        Hf=Sf*delta_x
        H2_2=altura_de_energia_1-Hf
        delta_H2=abs(H2_1-H2_2)
        delta_H2_final=delta_H2
        tirante_final_con_error += error
    return tirante_final_con_error

def preguntar_si_se_desea_hacer_mas_calculos_de_tirantes() :
    volver_a_calcular_tirantes = int(input("""Deseas calcular mas tirantes?
0.NO (si eliges esta opción el programa finalizará)
1.SI
Elige una opción: """))
    if volver_a_calcular_tirantes == 1 :
        if seccion_transversal == 1 :
            calcular_tirante_sección_rectangular()
        elif seccion_transversal == 2 :
            calcular_tirante_sección_trapecial()
    elif volver_a_calcular_tirantes == 0 :
        print('Gracias por usar mi Programa, atte: Jared DL')
    else :
        print('Seleciona una opcion válida!!!')
        preguntar_si_se_desea_hacer_mas_calculos_de_tirantes()

# FUNCIONES BÁSICAS:
# solera es lo mismo que base

def calcular_area_rectangulo(base,altura):
    area=base*altura
    return area

def calcular_perimetro_mojado_rectangulo(base,altura) :
    perimetro=2*altura+base
    return perimetro

def calcular_energia(gravedad,tirante,caudal,area) :
    energia=tirante+(caudal/area)**2/(2*gravedad)
    return energia

def calcular_pendiente(caudal,coef_manning,perimetro,area) :
    pendiente=((caudal*coef_manning*perimetro**(2/3))/area**(5/3))**2
    return pendiente

def calcular_area_trapecio(solera,tirante,Z_1):
    area=(solera+Z_1*tirante)*tirante
    return area

def calcular_perimetro_mojado_trapecio(solera,tirante,Z_1) :
    perimetro=2*tirante*(Z_1**2+1)**(1/2)+solera
    return perimetro

def calcular_altura_de_energía(altura_del_canal,tirante,caudal,area) :
    altura=altura_del_canal+tirante+(caudal/area)**2/(2*g_gravedad)
    return altura


inicio()