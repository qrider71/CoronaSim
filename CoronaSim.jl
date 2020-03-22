using CSV
using LsqFit
using Dates


function calc(dates::Array{Date}, data::Array{Float64}, n::Int64)

    function j_m_exp(t,p)
        J = Array{Float64}(undef, length(t), length(p))
        J[:,1] = exp.(p[2] .* t)       #df/dp[1]
        J[:,2] = t .* p[1] .* J[:,1]   #df/dp[2]
        J
    end

    function j_m_pow2(t,p)
        J = Array{Float64}(undef, length(t), length(p))
        J[:,1] = 2.0 .^ (t ./ p[2])       #df/dp[1]

        c = -log(2.0)*p[1]/(p[2] * p[2])
        J[:,2] =  t .* c .* J[:,1]        #df/dp[2]
        J
    end

    datesWindow = dates[end-n+1 : end]
    date1 = dates[end-n+1]

    x = map(d -> (d-date1).value, datesWindow)
    y = data[end-n+1 : end]

    m_exp(t, p) = p[1] * exp.(p[2] * t)
    p0_exp = [x[1], 0.5]
    fit_exp = curve_fit(m_exp, j_m_exp, x, y, p0_exp)

    m_pow2(t, p) = p[1] * 2.0 .^ (t/p[2])
    p0_pow2 = [x[1], 0.5]
    fit_pow2 = curve_fit(m_pow2, j_m_pow2, x, y, p0_pow2)

    confidence_intervals = confidence_interval(fit_pow2, 0.1)

    p = fit_pow2.param
    
    # y_fitted = m_pow2(x,p)
    # diff = 100.0*(abs.(y.-y_fitted))./y
    # println(diff)

    print(p)
    print("\t")
    println(confidence_intervals)
    println("---")

end

df = CSV.File("CoronaDataExport.csv")

dates = df.Date

dataDE = convert(Array{Float64}, df.Germany)
dataCH = convert(Array{Float64}, df.Switzerland)
dataIT = convert(Array{Float64}, df.Italy)
dataKR = convert(Array{Float64}, df.Korea)
dataSG = convert(Array{Float64}, df.Singapore)
dataUS = convert(Array{Float64}, df.USA)


print("Germany: ")
# calc(dates[1:end], dataDE, 10)
calc(dates[1:end], dataDE, 20)
# calc(dates[1:end], dataDE, 5)

print("Switzerland: ")
# calc(dates[1:end], dataCH, 15)
calc(dates[1:end], dataCH, 10)
# calc(dates[1:end], dataCH, 5)

print("Italy: ")
# calc(dates[1:end], dataIT, 15)
calc(dates[1:end], dataIT, 10)
# calc(dates[1:end], dataIT, 5)

print("South Korea: ")
# calc(dates[1:end], dataKR, 15)
calc(dates[1:end], dataKR, 10)
# calc(dates[1:end], dataKR, 5)

print("Singapore: ")
# calc(dates[1:end], dataSG, 15)
calc(dates[1:end], dataSG, 10)
# calc(dates[1:end], dataSG, 5)

print("USA: ")
# calc(dates[1:end], dataSG, 15)
calc(dates[1:end], dataUS, 10)
# calc(dates[1:end], dataSG, 5)


