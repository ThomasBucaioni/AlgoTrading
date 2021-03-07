#! /usr/bin/python
#----------
# This unusual and intriguing algorithm was originally invented 
# by Michael V. Klibanov, Professor, Department of Mathematics and Statistics,
# University of North  Carolina at Charlotte. It is published in the following
# paper:
# M.V. Klibanov, A.V. Kuzhuget and K.V. Golubnichiy, 
# "An ill-posed problem for the Black-Scholes equation 
# for a profitable forecast of prices of stock options on real market data", 
# Inverse Problems, 32 (2016) 015010.
#----------
# Script assumes it's called by crontab, at the opening of the market
#-----
import numpy as np
import pause, datetime
from bs4 import BeautifulSoup
import requests
# Quadratic interpolation of the bid and ask option prices, and linear interpolation in between (https://people.math.sc.edu/kellerlv/Quadratic_Interpolation.pdf)
def funcQuadraticInterpolationCoef(values): # There is 'scipy.interpolate.interp1d' too
    y = np.array(values)
    A = np.array([[1,0,0],[1,-1,1],[1,-2,4]])
    return np.linalg.solve(A,y) # https://en.wikipedia.org/wiki/Polynomial_regression
def funcUab(t,coef):
    return coef[2]*t**2 + coef[1]*t + coef[0]
def funcF(s, sa, sb, ua, ub):
    return (s-sb)*(ua-ub)/(sa-sb) + ub
# Initialize the volatility and option lists of 3 values
optionBid = [0] # dummy value to pop in the loop
optionAsk = [0] # dummy value to pop in the loop
volatility = [0] # dummy value to pop in the loop
# Initalization for the loop
Nt = 4 # even number greater than 2: 4, 6, ...
Ns = 2 # even number greater than 0: 2, 4, ...
twotau = 2 # not a parameter...
alpha = 0.01 # not a parameter...
dt = twotau / Nt # time grid step
dimA = ( (Nt+1)*(Ns+1), (Nt+1)*(Ns+1) ) # Matrix A dimensions
dimb = ( (Nt+1)*(Ns+1), 1 ) # Vector b dimensions
A = np.zeros( dimA ) # Matrix A
b = np.zeros( dimb ) # Vector b
portfolio = 1000000 # Money 'available'
securityMargin = 0.00083 # EMPIRICAL: needs to be adjusted when taking into account the transaction fees (should rise, see the article p.8)
# Wait 10mn after the opening of the market
datet = datetime.datetime.now()
datet = datetime.datetime(datet.year, datet.month, datet.day, datet.hour, datet.minute + 10)
pause.until(datet)
# Record the stock and option values and wait 10mn more
def funcRetrieveStockOptionVolatility():
    # Stock
    stock_data_url = "https://finance.yahoo.com/quote/MSFT?p=MSFT"
    stock_data_html = requests.get(data_url).content
    stock_content = BeautifulSoup(stock_data_html, "html.parser")
    stock_bid = content.find("td", {'class': 'Ta(end) Fw(600) Lh(14px)', 'data-test': "BID-value"})
    print(stock_bid)
    stock_ask = content.find("td", {'class': 'Ta(end) Fw(600) Lh(14px)', 'data-test': "ASK-value"})
    print(stock_ask)
    stockOptVol[0] = stock_bid.text.split()[0]
    stockOptVol[1] = stock_ask.text.split()[0]
    # Option
    option_data_url = "https://finance.yahoo.com/quote/MSFT/options?p=MSFT&date=1631836800"
    option_data_html = requests.get(option_data_url).content
    option_content = BeautifulSoup(option_data_html, "html.parser")
    call_option_table = content.find("table", {'class': 'calls W(100%) Pos(r) Bd(0) Pt(0) list-options'})
    calls = call_option_table.find_all("tr")[1:]
    it = 0
    for call_option in calls:
        it+=1
        print("it = ", it)
        if "in-the-money " in str(call_option):
            itm_calls.append(call_option)
            print("in the money")
            itm_put_data = []
            for td in BeautifulSoup(str(itm_calls[-1]), "html.parser").find_all("td"):
                itm_put_data.append(td.text)
            print(itm_put_data)
            if itm_put_data[0] == 'MSFT210917C00220000': # One single option
                stockOptVol[2] = float(itm_put_data[4])
                stockOptVol[3] = float(itm_put_data[5])
                stockOptVol[4] = float(itm_put_data[-1].strip('%'))
        else:
            otm_calls.append(call_option)
            print("out the money")
            print("bid = ", option_bid, "\nask = ", option_ask, "\nvol = ",option_vol)
    return stockOptVol
# Record option and volatility
stockOptVol = funcRetrieveStockOptionVolatility()
optionBid.append(stockOptVol[2])
optionAsk.append(stockOptVol[3])
optionVol.append(stockOptVol[4])
# Wait another 10mn to record a second value for the quadratic interpolation
datet = datetime.datetime.now()
datet = datetime.datetime(datet.year, datet.month, datet.day, datet.hour, datet.minute + 10)
pause.until(datet)
stockOptVol = funcRetrieveStockOptionVolatility()
optionBid.append(stockOptVol[2])
optionAsk.append(stockOptVol[3])
optionVol.append(stockOptVol[4])
tradeAtTimeTau = False
tradeAtTimeTwoTau = False
# Run the loop until 30mn before closure
datet = datetime.datetime.now()
datetend = datetime.datetime(datet.year, datet.month, datet.day, datet.hour + 6, datet.minute + 10)
while datet <= datetend:
    datet = datetime.datetime(datet.year, datet.month, datet.day, datet.hour, datet.minute + 10)
    optionBid.pop(0)
    optionAsk.pop(0)
    optionVol.pop(0)
    stockOptVol = funcRetrieveStockOptionVolatility()
    stockBid = stockOptVol[0]
    stockAsk = stockOptVol[1]
    optionBid.append(stockOptVol[2])
    optionAsk.append(stockOptVol[3])
    optionVol.append(stockOptVol[4])
    # Trade if required
    if tradeAtTimeTau == True or tradeAtTimeTwoTau == True: # sell
        if tradeAtTimeTau == True:
            portfolio += min(optionAsk[2],sellingPriceAtTimeTau) * 140 # sell 140 options bought 10mn ago
            tradeAtTimeTau = tradeAtTimeTwoTau
            sellingPriceAtTimeTau = sellingPriceAtTimeTwoTau
            sellingPriceAtTimeTwoTau = false
    else: # forecast the option when no trading
        # Interpolation
        coefa = funcQuadraticInterpolationCoef(optionAsk) # quadratic interpolation of the option ask price
        coefb = funcQuadraticInterpolationCoef(optionBid) # quadratic interpolation of the option bid price
        coefs = funcQuadraticInterpolationCoef(optionVol) # quadratic interpolation of the volatility sigma
        sa = stockAsk # stock ask price
        sb = stockBid # stock bid price
        ds = (sa - sb) / Ns # stock grid step
        for k in range (0, Ns+1): # fill the matrix and the vector
            for j in range (0, Nt+1):
                Atemp = np.zeros( dimA )
                btemp = np.zeros( dimb )
                print("k = {k}, j = {j}".format(k=k,j=j))
                if k == 0:
                    Atemp[ k*(Nt+1)+j, k*(Nt+1)+j ] = 1
                    btemp[ k*(Nt+1)+j ] = funcUab(j*dt,coefb)
                elif k == Ns:
                    Atemp[ k*(Nt+1)+j, k*(Nt+1)+j ] = 1
                    btemp[ k*(Nt+1)+j ] = funcUab(j*dt,coefa)
                elif j == 0:
                    Atemp[ k*(Nt+1)+j, k*(Nt+1)+j ] = 1
                    btemp[ k*(Nt+1)+j ] = funcF( k*ds+sb, sa, sb, funcUab(j*dt,coefa), funcUab(j*dt,coefb) )
                elif j == Nt: # do nothing
                    pass
                else: # main case
                    akj = 0.5*(255*13*3)* funcUab(j*dt, coefs)**2 * (k*ds + sb)**2
                    dts = (twotau-dt)/Nt * (sa-sb-ds)/Ns
                    #----------
                    #----- Integral of the generator L
                    #----------
                    #----- time derivative
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j+1) ] =   dts / dt**2 # k,j+1 ~ k,j+1
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j-1) ] =   dts / dt**2 # k,j-1 ~ k,j-1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j-1) ] = - dts / dt**2 # k,j+1 ~ k,j-1
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j+1) ] = - dts / dt**2 # k,j-1 ~ k,j+1
                    #----------
                    #----- stock derivative
                    #----------
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] =      akj**2 * dts / ds**4 # k+1,j ~ k+1,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] =  4 * akj**2 * dts / ds**4 # k,j ~ k,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] =      akj**2 * dts / ds**4 # k-1,j ~ k-1,j
                    #-----
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] = -2 * akj**2 * dts / ds**4 # k+1,j ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] = -2 * akj**2 * dts / ds**4 # k,j ~ k+1,j
                    #-----
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] = -2 * akj**2 * dts / ds**4 # k-1,j ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] = -2 * akj**2 * dts / ds**4 # k,j ~ k-1,j
                    #-----
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] =      akj**2 * dts / ds**4 # k+1,j ~ k-1,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] =      akj**2 * dts / ds**4 # k-1,j ~ k+1,j
                    #----------
                    #----- time and stock derivatives
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+1)*(Nt+1)+(j+0) ] =      akj    * dts / (dt*ds**2) # k,j+1 ~ k+1,j
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+1) ] =      akj    * dts / (dt*ds**2) # k+1,j ~ k,j+1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j+0) ] =    - akj    * dts / (dt*ds**2) # k,j-1 ~ k+1,j
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j-1) ] =    - akj    * dts / (dt*ds**2) # k+1,j ~ k,j-1
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j+0) ] = -2 * akj    * dts / (dt*ds**2) # k,j+1 ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+1) ] = -2 * akj    * dts / (dt*ds**2) # k,j ~ k,j+1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j+0) ] =  2 * akj    * dts / (dt*ds**2) # k,j-1 ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j-1) ] =  2 * akj    * dts / (dt*ds**2) # k,j ~ k,j-1
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k-1)*(Nt+1)+(j+0) ] =      akj    * dts / (dt*ds**2) # k,j+1 ~ k-1,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+1) ] =      akj    * dts / (dt*ds**2) # k-1,j ~ k,j+1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j+0) ] =    - akj    * dts / (dt*ds**2) # k,j-1 ~ k-1,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j-1) ] =    - akj    * dts / (dt*ds**2) # k-1,j ~ k,j-1
                    #----------
                    #----------
                    #----- Regularisation term - using alpha = 0.01
                    #----------
                    #---------- 
                    #----- H2 norm: 0 derivative
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] += alpha # k,j ~ k,j
                    #-----
                    coef = funcF( k*ds+sb, sa, sb, funcUab(j*dt,coefa), funcUab(j*dt,coefb) )
                    btemp[ (k+0)*(Nt+1)+(j+0) ] +=  alpha * 2 * coef
                    #----------
                    #----- H2 norm: time derivative
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j+1) ] +=  alpha / dt**2 # k,j+1 ~ k,j+1
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j-1) ] +=  alpha / dt**2 # k,j-1 ~ k,j-1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j-1) ] += -alpha / dt**2 # k,j+1 ~ k,j-1
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j+1) ] += -alpha / dt**2 # k,j-1 ~ k,j+1
                    #-----
                    coef = ( funcF( k*ds+sb, sa, sb, funcUab((j+1)*dt,coefa), funcUab((j+1)*dt,coefb) ) \
                    -        funcF( k*ds+sb, sa, sb, funcUab((j-1)*dt,coefa), funcUab((j-1)*dt,coefb) ) ) / dt
                    btemp[ (k+0)*(Nt+1)+(j+1) ] +=   alpha * 2 * coef
                    btemp[ (k+0)*(Nt+1)+(j-1) ] += - alpha * 2 * coef
                    #----------
                    #----- H2 norm: stock derivative
                    #----------
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] +=  alpha / ds**2 # k+1,j ~ k+1,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] +=  alpha / ds**2 # k-1,j ~ k-1,j
                    #-----
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] += -alpha / ds**2 # k+1,j ~ k-1,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] += -alpha / ds**2 # k-1,j ~ k+1,j
                    #-----
                    coef = ( funcUab(j*dt,coefa) - funcUab(j*dt,coefb) ) / (sa - sb)
                    btemp[ (k+1)*(Nt+1)+(j+0) ] +=   alpha * 2 * coef
                    btemp[ (k-1)*(Nt+1)+(j+0) ] += - alpha * 2 * coef
                    #----------
                    #----- H2 norm: stock and time derivative
                    #----------
                    Atemp[ (k+1)*(Nt+1)+(j+1), (k+1)*(Nt+1)+(j+1) ] +=  alpha / (ds*dt) # k+1,j+1 ~ k+1,j+1
                    Atemp[ (k-1)*(Nt+1)+(j+1), (k-1)*(Nt+1)+(j+1) ] +=  alpha / (ds*dt) # k-1,j+1 ~ k-1,j+1
                    Atemp[ (k-1)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j-1) ] +=  alpha / (ds*dt) # k-1,j-1 ~ k-1,j-1
                    Atemp[ (k+1)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j-1) ] +=  alpha / (ds*dt) # k+1,j-1 ~ k+1,j-1
                    #----------
                    Atemp[ (k+1)*(Nt+1)+(j+1), (k-1)*(Nt+1)+(j+1) ] += -alpha / (ds*dt) # k+1,j+1 ~ k-1,j+1
                    Atemp[ (k+1)*(Nt+1)+(j+1), (k+1)*(Nt+1)+(j-1) ] += -alpha / (ds*dt) # k+1,j+1 ~ k+1,j-1
                    Atemp[ (k+1)*(Nt+1)+(j+1), (k-1)*(Nt+1)+(j-1) ] +=  alpha / (ds*dt) # k+1,j+1 ~ k-1,j-1
                    #-----
                    Atemp[ (k-1)*(Nt+1)+(j+1), (k+1)*(Nt+1)+(j+1) ] += -alpha / (ds*dt) # k-1,j+1 ~ k+1,j+1
                    Atemp[ (k+1)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j+1) ] += -alpha / (ds*dt) # k+1,j-1 ~ k+1,j+1
                    Atemp[ (k-1)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j+1) ] +=  alpha / (ds*dt) # k-1,j-1 ~ k+1,j+1
                    #----------
                    Atemp[ (k-1)*(Nt+1)+(j+1), (k+1)*(Nt+1)+(j-1) ] +=  alpha / (ds*dt) # k-1,j+1 ~ k+1,j-1
                    Atemp[ (k-1)*(Nt+1)+(j+1), (k-1)*(Nt+1)+(j-1) ] += -alpha / (ds*dt) # k-1,j+1 ~ k-1,j-1
                    #-----
                    Atemp[ (k+1)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j+1) ] +=  alpha / (ds*dt) # k+1,j-1 ~ k-1,j+1
                    Atemp[ (k-1)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j+1) ] += -alpha / (ds*dt) # k-1,j-1 ~ k-1,j+1
                    #----------
                    Atemp[ (k+1)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j-1) ] += -alpha / (ds*dt) # k+1,j-1 ~ k-1,j-1
                    #-----
                    Atemp[ (k-1)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j-1) ] += -alpha / (ds*dt) # k-1,j-1 ~ k+1,j-1
                    #----------
                    coef = ( funcUab((j+1)*dt,coefa) - funcUab((j+1)*dt,coefb) \
                    -        funcUab((j-1)*dt,coefa) + funcUab((j-1)*dt,coefb) ) / (dt * (sa - sb))
                    btemp[ (k+1)*(Nt+1)+(j+1) ] +=   alpha * 2 * coef / (ds*dt)
                    btemp[ (k-1)*(Nt+1)+(j+1) ] += - alpha * 2 * coef / (ds*dt)
                    btemp[ (k-1)*(Nt+1)+(j-1) ] += - alpha * 2 * coef / (ds*dt)
                    btemp[ (k+1)*(Nt+1)+(j-1) ] +=   alpha * 2 * coef / (ds*dt)
                    #----------
                    #----- H2 norm: stock second derivative
                    #----------
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j+1) ] +=      alpha / dt**4 # k,j+1 ~ k,j+1
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] +=  4 * alpha / dt**4 # k,j ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j-1) ] +=      alpha / dt**4 # k,j-1 ~ k,j-1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j+0) ] += -2 * alpha / dt**4 # k,j+1 ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+1) ] += -2 * alpha / dt**4 # k,j ~ k,j+1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j+1), (k+0)*(Nt+1)+(j-1) ] +=      alpha / dt**4 # k,j+1 ~ k,j-1
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j+1) ] +=      alpha / dt**4 # k,j-1 ~ k,j+1
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j-1) ] += -2 * alpha / dt**4 # k,j ~ k,j-1
                    Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j+0) ] += -2 * alpha / dt**4 # k,j-1 ~ k,j
                    #----------
                    #----- H2 norm: time second derivative
                    #----------
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] +=      alpha / ds**4 # k+1,j ~ k+1,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] +=  4 * alpha / ds**4 # k,j ~ k,j
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] +=      alpha / ds**4 # k-1,j ~ k-1,j
                    #-----
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] += -2 * alpha / ds**4 # k+1,j ~ k,j
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] += -2 * alpha / ds**4 # k,j ~ k+1,j
                    #-----
                    Atemp[ (k+1)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] +=      alpha / ds**4 # k,j ~ k,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] +=      alpha / ds**4 # k,j ~ k,j
                    #-----
                    Atemp[ (k+0)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] += -2 * alpha / ds**4 # k,j ~ k-1,j
                    Atemp[ (k-1)*(Nt+1)+(j+0), (k+0)*(Nt+1)+(j+0) ] += -2 * alpha / ds**4 # k-1,j ~ k,j
                    #----------
                    coef = ( funcF( k*ds+sb, sa, sb, funcUab((j+1)*dt,coefa), funcUab((j+1)*dt,coefb) ) \
                    -    2 * funcF( k*ds+sb, sa, sb, funcUab((j+0)*dt,coefa), funcUab((j+0)*dt,coefb) ) \
                    +        funcF( k*ds+sb, sa, sb, funcUab((j-1)*dt,coefa), funcUab((j-1)*dt,coefb) ) ) / dt**2
                    btemp[ (k+0)*(Nt+1)+(j+1) ] +=   alpha * 2 * coef / dt**2
                    btemp[ (k+0)*(Nt+1)+(j+0) ] += - alpha * 4 * coef / dt**2
                    btemp[ (k+0)*(Nt+1)+(j-1) ] +=   alpha * 2 * coef / dt**2
                    #----------
                    #----------
                    #----- Boundary de-computation
                    #----------
                    if k+1 == Ns:
                        Atemp[ (k+1)*(Nt+1)+(j+0), (k+1)*(Nt+1)+(j+0) ] = 0 # k+1,j ~ k+1,j
                        Atemp[ (k+1)*(Nt+1)+(j+1), (k+1)*(Nt+1)+(j+1) ] = 0 # k+1,j+1 ~ k+1,j+1
                        Atemp[ (k+1)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j-1) ] = 0 # k+1,j-1 ~ k+1,j-1
                        btemp[ (k+1)*(Nt+1)+(j+0) ] = 0 # k+1,j
                        btemp[ (k+1)*(Nt+1)+(j+1) ] = 0 # k+1,j+1
                        btemp[ (k+1)*(Nt+1)+(j-1) ] = 0 # k+1,j-1
                    if k-1 == 0:
                        Atemp[ (k-1)*(Nt+1)+(j+0), (k-1)*(Nt+1)+(j+0) ] = 0 # k-1,j ~ k-1,j
                        Atemp[ (k-1)*(Nt+1)+(j+1), (k-1)*(Nt+1)+(j+1) ] = 0 # k-1,j+1 ~ k-1,j+1
                        Atemp[ (k-1)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j-1) ] = 0 # k-1,j-1 ~ k-1,j-1
                        btemp[ (k-1)*(Nt+1)+(j+0) ] = 0 # k-1,j
                        btemp[ (k-1)*(Nt+1)+(j+1) ] = 0 # k-1,j+1
                        btemp[ (k-1)*(Nt+1)+(j-1) ] = 0 # k-1,j-1
                    if j-1 == 0:
                        Atemp[ (k+0)*(Nt+1)+(j-1), (k+0)*(Nt+1)+(j-1) ] = 0 # k,j-1 ~ k,j-1
                        Atemp[ (k+1)*(Nt+1)+(j-1), (k+1)*(Nt+1)+(j-1) ] = 0 # k+1,j-1 ~ k+1,j-1
                        Atemp[ (k-1)*(Nt+1)+(j-1), (k-1)*(Nt+1)+(j-1) ] = 0 # k-1,j-1 ~ k-1,j-1
                        btemp[ (k+0)*(Nt+1)+(j-1) ] = 0 # k,j-1
                        btemp[ (k+1)*(Nt+1)+(j-1) ] = 0 # k+1,j-1
                        btemp[ (k-1)*(Nt+1)+(j-1) ] = 0 # k-1,j-1
                    #----------
                    pass
                print("-----")
                print("Atemp = ")
                print(Atemp)
                print("-----")
                print("btemp = ")
                print(btemp)
                print("-----")
                print("-----")
                A = A + Atemp
                b = b + btemp
                print("-----")
                print("A = ")
                print(A)
                print("-----")
                print("b = ")
                print(b)
                print("-----")
                print("-----")
                input("Press Enter to continue...")
        # Conjugate gradient algorithm: https://en.wikipedia.org/wiki/Conjugate_gradient_method
        x = np.zeros(N).reshape(N,1)
        r = b - np.matmul(A,x)
        p = r
        rsold = np.dot(r.transpose(),r)
        for i in range(len(b)):
            Ap = np.matmul(A,p)
            alpha = rsold / np.matmul(p.transpose(),Ap)
            x = x + alpha * p
            r = r - alpha * Ap
            rsnew = np.dot(r.transpose(),r)
            if np.sqrt(rsnew) < 1e-16:
                break
            p = r + (rsnew / rsold) * p
            rsold = rsnew
            print("it = ", i)
            print("rsold = ", rsold)
        # Trading strategy
        sm = (sa + sb)/2
        if x[Ns/2*(Nt+1)+Nt/2] >= optionAsk[0] + securityMargin:
            tradeAtTimeTau = True
            sellingPriceAtTimeTau = x[Ns/2*(Nt+1)+Nt/2]
            portfolio -= 140 * optionAsk # buy 140 options
        if x[Ns/2*(Nt+1)+Nt] >= optionAsk[0] + securityMargin:
            tradeAtTimeTwoTau = True
            sellingPriceAtTimeTwoTau = x[Ns/2*(Nt+1)+Nt]
            portfolio -= 140 * optionAsk # buy 140 options
            pause.until(datet)
    # Wait 10mn before the next loop
    pause.until(datet)
    datet = datetime.datetime.now()
# Time should be around 20mn before closure
datet = datetime.datetime(datet.year, datet.month, datet.day, datet.hour, datet.minute + 10)
if tradeAtTimeTau == True: # sell
    stockOptVol = funcRetrieveStockOptionVolatility()
    optionAsk.pop(0)
    optionAsk.append(stockOptVol[3])
    portfolio += min(optionAsk[2],sellingPriceAtTimeTau) * 140
# Wait 10mn more to sell the last options
pause.until(datet) # it should be around 10mn before closure
if tradeAtTimeTwoTau == True: # sell
    stockOptVol = funcRetrieveStockOptionVolatility()
    optionAsk.pop(0)
    optionAsk.append(stockOptVol[3])
    portfolio += min(optionAsk[2],sellingPriceAtTimeTwoTau) * 140
# Market closure
