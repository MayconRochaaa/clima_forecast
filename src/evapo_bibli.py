import pandas as pd
import numpy as np
import urllib3
import requests
import urllib
import time

from pandas import Timestamp

# Pacote parcial das funções meteriologicas do pacote simET (escrito em R).
#Fonte: https://rdrr.io/cran/simET/

# ========================================================================================================================
def cal_atmosphericPressure(elevation):
    """Calcula a pressão atmosférica com base na elevação

    Args:
        elevation (float): elevação a nivel do mar (metros)

    Returns:
        float: Pressão atmosférica (kPa)
    """
        
    P = 101.3*((293-0.0065*elevation)/293)**5.26
    
    return P

# ========================================================================================================================
def cal_slopeOfSaturationVapourPressureCurve(Tem):
    """Cálculo da inclinação da curva de pressão de vapor de saturação

    Args:
        Tem (float): Temperatura (ºC)

    Returns:
        float: inclinação da curva de pressão de vapor de saturação
    """

    slopeOfSaturationVapourPressureCurve = (4098*(0.6108*np.exp((17.27*Tem)/(Tem+237.3))))/(Tem+237.3)**2

    return   slopeOfSaturationVapourPressureCurve

# ========================================================================================================================
def cal_psychrometriCconstant(atmospheric_pressure):
    """Cálculo da constante psicrométrica

    Args:
        atmospheric_pressure (float): pressão atmosférica (kPa)

    Returns:
        float: constante psicrométrica (kPa/ºC)
    """
        
    gamma = 0.665**10**(-3)*atmospheric_pressure
    
    return gamma

# ========================================================================================================================
def cal_saturationVapourPressure(Tem):
    """Calculo da pressão de vapor de saturação

    Args:
        Tem (float): Temperatura (ºC)

    Returns:
        float: pressão de vapor de saturação (kPa)
    """
        
    saturation_vapour_pressure = 0.6108*np.exp((17.27*Tem)/(Tem+237.3))

    return saturation_vapour_pressure

# ========================================================================================================================
def cal_meanSaturationVapourPressure(Tmax, Tmin):
    """Calculo a pressão média de vapor de saturação

    Args:
        Tmax (float): Temperatura máxima do dia (ºC)
        Tmin (float): Temperatura miníma do dia (ºC)

    Returns:
        float: pressão média de vapor de saturação (kPa)
    """
    
    mean_saturation_vapour_pressure = (cal_saturationVapourPressure(Tmax)+cal_saturationVapourPressure(Tmin))/2
    
    return mean_saturation_vapour_pressure

# ========================================================================================================================
def cal_ActualVapourPressure_from_RHmean(RHmean,Tmax,Tmin):
    """Cálculo da pressão de vapor real derivada de RHmean

    Args:
        RHmean (float): umidade relativa média
        Tmax (float): temperatura máxima diária
        Tmin (float): temperatura mínima diária
        
    Returns:
        float: pressão de vapor real (kPa)
    """
    
    ActualVapourPressure = (RHmean/100)*((cal_saturationVapourPressure(Tmax)+cal_saturationVapourPressure(Tmin))/2)
    
    return ActualVapourPressure

# ========================================================================================================================
def convert_angert_to_radian(anger):
    """Converte angulo de graus para radianos

    Args:
        anger (float): Angulo em graus

    Returns:
        float: Angulo em radianos
    """
    radian = float(anger)*np.pi/180
    
    return radian

# ========================================================================================================================
def get_day_of_year(x):
    """Obtém o dia do ano de uma data x fornecida

    Args:
        x (timestamp): Data no formato YYYY/mm/dd

    Raises:
        ValueError1: Dia fora do intervalo contido em um ano
        ValueError2: Unidade fora da do esperado

    Returns:
        int: dia do ano
    """
    
    if isinstance(x, list) and all(isinstance(timepoint, Timestamp) for timepoint in x):
        doy = [int(timepoint.strftime("%j")) for timepoint in x]
    elif isinstance(x, Timestamp):
        doy = int(x.strftime("%j"))
    elif isinstance(x, (int, float)):
        if x <= 0 or x > 366:
            raise ValueError("Day of the year is not between 1-366")
        doy = x
    else:
        raise ValueError("Invalid input type. Please provide a Timestamp, a list of Timestamps, or a numeric value.")
    
    return doy

# ========================================================================================================================
def cal_inverseRelativeDistance_Earth_sun(J):
    """Calculo do fator de correção  da distância relativa Terra-Sol

    Args:
        J (timestamp): data 

    Returns:
        float: fator de correção para a distância relativa entre a Terra e o Sol
    """
    
    doy = get_day_of_year(J)
    
    inverseRelativeDistance_Earth_sun = 1+0.033*np.cos((2*np.pi/365)*doy)
    
    return inverseRelativeDistance_Earth_sun


# ========================================================================================================================
def cal_solarDeclination_in_FAO(J):
    """Cálculo da declinação solar com o método FAO56

    Args:
        J (timestamp): data

    Returns:
        float: fator de declinação solar
    """
    
    doy = get_day_of_year(J)
    
    solarDeclination = 0.409*np.sin((2*np.pi/365)*doy-1.39)
    
    return solarDeclination

# ========================================================================================================================
def cal_sunsetHourAngle(lat, solar_declination):
    """Calculo do ângulo da hora do pôr do sol


    Args:
        lat (float): latitude
        solar_declination (float): fator de declinação solar

    Returns:
        float: ângulo da hora do pôr do sol (radiano)
    """
    
    sunsetHourAngle = np.arccos(-np.tan(lat)*np.tan(solar_declination))

    return sunsetHourAngle
    
# ========================================================================================================================
def cal_extraterrestrialRadiation_for_daily(J,lat):
    """A radiação extraterrestre, Ra, para cada dia do ano e para diferentes latitudes pode ser estimada a partir da constante solar, da declinação solar e da época do ano.

    Args:
        J (timestamp): data
        lat (float): latitude

    Returns:
        float: radiação extraterrestre diária (MJ m-2 dia-1)
    """
    
    Gsc = 0.082
    dr = cal_inverseRelativeDistance_Earth_sun(J)
    solar_declination = cal_solarDeclination_in_FAO(J)
    ws = cal_sunsetHourAngle(lat,solar_declination)
    Ra = 24*60*0.082*dr*(ws*np.sin(lat)*np.sin(solar_declination)+np.cos(lat)*np.cos(solar_declination)*np.sin(ws))/np.pi
    
    return Ra

# ========================================================================================================================
def cal_skySolarRadiation_withas_elevation(z, Ra):
    """Cálculo da radiação de céu claro, Rso.

    Args:
        z (float): elevação da estação acima do nível do mar [m].
        Ra (float): Radiação extraterrestre [MJ m-2 dia-1]

    Returns:
        float: radiação solar de céu claro [MJ m-2 dia-1]
    """
    Rso = (0.75+2*10**(-5)*z)*Ra
    
    return Rso

# ========================================================================================================================
def cal_netSolarRadiation(alpha, Rs):
    """A radiação líquida de ondas curtas resultante do equilíbrio entre a radiação solar recebida e refletida.

    Args:
        alpha (float): coeficiente de reflexão do dossel
        Rs (float): radiação solar incidente [MJ m-2 dia-1]

    Returns:
        float:  radiação solar líquida [MJ m-2 dia-1]
    """
    
    Rns = (1-alpha)*Rs
    
    return Rns

# ========================================================================================================================
def convert_degreesCelsius_to_Kelvin(degrees_Celsius):
    """Conversão de graus Celsius para Kelvin

    Args:
        degrees_Celsius (float): temperatura (ºC)

    Returns:
        float: temperatura (Kelvin)
    """
    
    Kelvin = degrees_Celsius+273.16
    
    return Kelvin

# ========================================================================================================================
def cal_netLongwaveRadiation(TKmax,TKmin,ea,Rs,Rso):
    """Cálculo da radiação líquida de ondas longas Rnl

    Args:
        TKmax (float): temperatura absoluta máxima durante o período de 24 horas [K]
        TKmin (float): temperatura absoluta mínima durante o período de 24 horas [K]
        ea (float): pressão de vapor real [kPa].
        Rs (float): radiação solar medida ou calculada [MJ m-2 dia-1]
        Rso (float): radiação de céu claro calculada [MJ m-2 dia-1].

    Returns:
        float: radiação líquida de ondas longas de saída [MJ m-2 dia-1]
    """
    
    delta = 4.903*10**(-9)#Constante de Stefan-Boltzmann MJ K-4 m-2 day-1
    Rnl = delta*((TKmax**4+TKmin**4)/2)*(0.34-0.14*np.sqrt(ea))*(1.35*(Rs/Rso)-0.35)
    
    return Rnl

# ========================================================================================================================
def cal_netRadiation(Rns, Rnl):
    """Cálculo da radiação líquida Rn

    Args:
        Rns (float): radiação líquida de ondas curtas recebida [MJ m-2 dia-1]
        Rnl (float): adiação líquida de ondas longas de saída [MJ m-2 dia-1]

    Returns:
        float: radiação liquida [MJ m-2 dia-1]
    """
    Rn = Rns - Rnl
    
    return Rn

# ========================================================================================================================
def cal_ET0_from_PM(delta,Rn,G,gamma,Tem,u2,es,ea):
    """cálculo da evapotranspiração de referência a partir do método Penman-Monteith

    Args:
        delta (float): curva de pressão de vapor inclinada (kPa °C)
        Rn (float): Radiação líquida na superfície da cultura [MJ m-2 dia-1]
        G (float): densidade de fluxo de calor no solo [MJ m-2 dia-1].
        gamma (float): constante psicrométrica (kPa °C).
        Tem (float): temperatura do ar a 2 m de altura [°C].
        u2 (float): velocidade do vento a 2 m de altura [m s-1].
        es (float): pressão de vapor de saturação [kPa].
        ea (float): pressão de vapor real [kPa].
        
    Returns:
        float: evapotranspiração de referência [mm dia-1]
    """
    

    ET0 = (0.408*delta*(Rn-G)+gamma*(900/(Tem+273))*u2*(es-ea))/(delta+gamma*(1+0.34*u2))

    return ET0

# ========================================================================================================================
def cal_ET0_from_PM_for_daily(Latitude,Altitude,J,Tmax,Tmin,Rs,RHmean,Wind):
    """Cálculo da evapotranspiração de referência de Penman-Monteith por dia

    Args:
        Latitude (float): latitude
        Altitude (float): altitude a nivel do mar (m)
        J (timestamp): data
        Tmax (float): temperatura máxima diária do ar (ºC)
        Tmin (float): temperatura mínima diária do ar (ºC)
        Rs (float): Radiação solar [MJ m-2 d-1]
        RHmean (float): umidade relativa média diária
        Wind (float): velocidade do vento a 2 m de altura [m s-1].

    Returns:
        float: evapotranspiração de referência (mm/dia)
    """
    
    Tmean=(Tmax+Tmin)/2
    P=cal_atmosphericPressure(Altitude)
    Delta=cal_slopeOfSaturationVapourPressureCurve(Tmean)
    gamma=cal_psychrometriCconstant(P)
    es=cal_meanSaturationVapourPressure(Tmax ,Tmin)
    ea=cal_ActualVapourPressure_from_RHmean(RHmean,Tmax,Tmin)
    Deficit=es-ea
    # dr=cal_inverseRelativeDistance_Earth_sun(J)
    # Solar_D=cal_solarDeclination_in_FAO(J)
    Lat=convert_angert_to_radian(Latitude)
    # ws=cal_sunsetHourAngle(Lat,Solar_D)
    Ra=cal_extraterrestrialRadiation_for_daily(J,Lat)
    # Nmax=cal_daylightHours(ws)
    # Rs=cal_solarRadiation(as=0.25,bs=0.5,n=Na,N=Nmax,Ra=Ra)
    Rso=cal_skySolarRadiation_withas_elevation(z=Altitude,Ra=Ra)
    Rns=cal_netSolarRadiation(alpha=0.23,Rs=Rs)
    TKmax=convert_degreesCelsius_to_Kelvin(Tmax)
    TKmin=convert_degreesCelsius_to_Kelvin(Tmin)
    Rnl=cal_netLongwaveRadiation(TKmax,TKmin,ea,Rs,Rso)
    Rn=cal_netRadiation(Rns,Rnl)
    G=0
    ET0=cal_ET0_from_PM(Delta,Rn,G,gamma,Tmean,Wind,es,ea)
    
    return ET0

# ========================================================================================================================
def estimate_ET0_with_TmaxAndTmin(Tmean,Tmax,Tmin,Ra):
    """Estima a evapotranspiração com Tmax e Tmin

    Args:
        Tmean (_type_): temperatura média do dia (ºC)
        Tmax (_type_): temperatura maxima do dia (ºC)
        Tmin (_type_): temperatura minima do dia (ºC)
        Ra (_type_): radiação extraterrestre [mm dia-1]
        
    Returns:
        float: evapotranspiração de referência (mm dia-1).
    """
  
    ET0 = 0.0023*(Tmean+17.8)*(Tmax-Tmin)**0.5*Ra
    
    return ET0

# ========================================================================================================================

def estimate_Rs_from_airTemDiff(Ra, Tmax, Tmin, locations):
    """Estima dados de radiação solar derivados de diferenças de temperatura do ar


    Args:
        Ra (float): radiação extraterrestre [MJ m-2 d-1].
        Tmax (float): temperatura máxima do ar (ºC).
        Tmin (float): temperatura mínima do ar (ºC).
        locations (string): região de interior ou costeira

    Returns:
        float: radiação solar [MJ m-2 d-1]
    """
    
    # if locations == "interior":
    #     k_Rs = 0.16
    # elif locations == "coastal":
    #     k_Rs = 0.19
    # k_Rs = ifelse(locations=="interior",k_Rs<-0.16,
    #               ifelse(locations=="coastal",k_Rs<-0.19,
    #                      NA))
    k_Rs = 0.16 if locations == "interior" else (0.19 if locations == "coastal" else None)
    Rs = k_Rs * (Tmax - Tmin)**0.5 * Ra
    return Rs

# ========================================================================================================================
def make_remote_request(url: str, params: dict):
   """
   Makes the remote request
   Continues making attempts until it succeeds
   """

   count = 1
   while True:
       try:
           response = requests.get((url + urllib.parse.urlencode(params)))
       except (OSError, urllib3.exceptions.ProtocolError) as error:
           print('\n')
           print('*' * 20, 'Error Occured', '*' * 20)
           print(f'Number of tries: {count}')
           print(f'URL: {url}')
           print(error)
           print('\n')
           count += 1
           continue
       break

   return response

# ========================================================================================================================

def elevation_function(x):
    """Realiza a requisição de dados de elevação com base nos dados de latitude e longitude

    Args:
        x (DataFrame): DataFrame com dados de latitude e longitude (somente)
        
    Returns:
        pd.Series: Dados de elevação
    """
    
    url = 'https://api.opentopodata.org/v1/srtm30m?'
    params = {'locations': f"{x[0]},{x[1]}"}
    result = make_remote_request(url, params)
    time.sleep(1)
    return result.json()['results'][0]['elevation']

# ========================================================================================================================

def get_elevations_for_stations(df):
    """Calcula a elevação de cada estação a partir da API

    Args:
        df (DataFrame): DataFrame contendo dados de latitude e longitude e estação (somente)

    Returns:
        pd.Series: Dados de elevação do dataframe fornecido
    """
    unique_stations = df['estacao'].unique()
    elevations = []

    for station in unique_stations:
        station_data = df[df['estacao'] == station].iloc[0]
        coordinates = (station_data['lat'], station_data['lon'])
        elevation = elevation_function(coordinates)
        elevations.append(elevation)

    return pd.Series(elevations, index=unique_stations)