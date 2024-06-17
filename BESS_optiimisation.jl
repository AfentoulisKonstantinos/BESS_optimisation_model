using CSV, DataFrames, Dates, JuMP, HiGHS

# Define a custom date format
date_format = DateFormat("dd/mm/yy HH:MM")

# Load the data from the CSV file
HourlyData = CSV.read("half_hourly_data_for_all_markets.csv", DataFrame)
HourlyData.timestamp = [DateTime(row.timestamp, date_format) for row in eachrow(HourlyData)]

# Get unique dates
unique_dates = unique(Date.(HourlyData.timestamp))

# Load parameters from CSV
parameters = CSV.read("BESS_parameters.csv", DataFrame)

# Convert parameters DataFrame to a Dict
param_dict = Dict(row.parameter => row.value for row in eachrow(parameters))

# Define battery parameters from the parameter dictionary
## Set the maximum and the minimum charge and discharge rate, in MW
max_charge_rate = param_dict["max_charge_rate_MW"]
max_discharge_rate = param_dict["max_discharge_rate_MW"]
min_charge_rate = param_dict["min_charging_rate_MW"]
min_discharge_rate = param_dict["min_discharging_rate_MW"]

## Set the maximum and the minimum level of the stored energy in the battery - avoid damaging, in MWh
max_storage_capacity_limit = param_dict["max_SOE"] * param_dict["max_storage_capacity_MWh"]
min_storage_capacity_limit = param_dict["min_SOE"] * param_dict["max_storage_capacity_MWh"]

## Set the total capacity of the battery - avoid damaging, in MWh
max_storage_capacity = param_dict["max_storage_capacity_MWh"]

## Set the initial and the final energy stored in the battery in MWh
initial_storage = param_dict["initial_storage_level_pu"] * param_dict["max_storage_capacity_MWh"]
final_storage = param_dict["final_storage_level_pu"] * param_dict["max_storage_capacity_MWh"]

## Set the efficiency parameters
charging_efficiency = 1 - param_dict["charging_efficiency"]
discharging_efficiency = 1 - param_dict["discharging_efficiency"]

## Calculate the annual cycle limit - This is calculated as the total lifetime cycles divided by the total lifetime years
cycles_per_year = param_dict["lifetime_cycles"] / param_dict["lifetime_years"]

## Set the duration of the optimisation horizon
optimisation_horizon_years  = param_dict["optimisation_horizon"]

# Set the duration of each settlement period
settlement_period_duration = param_dict["settlement_period_duration"]

# Initialize cumulative charge, discharge power, and profits
global cumulative_charged_power = 0.0
global cumulative_discharged_power = 0.0
global cumulative_profit = 0.0

# Initialize a DataFrame to store the results
results_df = DataFrame()

# Create a model, use HiGHS as a solver
model = Model(HiGHS.Optimizer)

# Set the maximum time limit in seconds - float64 
set_optimizer_attribute(model, "time_limit", 100.0)

# Calculate the number of settlement periods
daily_settlement_periods = count(row -> Date(row.timestamp) == unique_dates[1], eachrow(HourlyData))

# Define the set of time periods and days
time_hours = 1:48 

print("time_hours", time_hours, "daily_settlement_periods : ", daily_settlement_periods)
days = 1:length(unique_dates)

# Decision variables

## Decision variables for Market 1
@variable(model, charging_power_m1[days, time_hours] >= 0)
@variable(model, discharging_power_m1[days, time_hours] >= 0)

## Decision variables for Market 2
@variable(model, charging_power_m2[days, time_hours] >= 0)
@variable(model, discharging_power_m2[days, time_hours] >= 0)

## Decision variables for Market 3
@variable(model, charging_power_m3[days,time_hours] >= 0)
@variable(model, discharging_power_m3[days,time_hours] >= 0)

## Decision variables for participation in all three markets
@variable(model, charge_decision[days, time_hours], Bin)
@variable(model, discharge_decision[days, time_hours], Bin)
@variable(model, charging_power[days, time_hours] >= 0)
@variable(model, discharging_power[days, time_hours] >= 0)

## Decision variable for the storage level
@variable(model, storage_level[days, time_hours] >= 0)

## Variable to track cumulative energy throughput
@variable(model, cumulative_energy_throughput >= 0)

# Objective function
@objective(model, Max,
   sum(discharging_power_m1[d, t] * discharging_efficiency * HourlyData[(d-1)*48 + t, :Market1Price] for d in days, t in time_hours) -
   sum(charging_power_m1[d, t] / charging_efficiency * HourlyData[(d-1)*48 + t, :Market1Price] for d in days, t in time_hours) +
   sum(discharging_power_m2[d, t] * discharging_efficiency* HourlyData[(d-1)*48 + t, :Market2Price] for d in days, t in time_hours) -
   sum(charging_power_m2[d, t] / charging_efficiency * HourlyData[(d-1)*48 + t, :Market2Price] for d in days, t in time_hours) +
   sum(discharging_power_m3[d,t] * discharging_efficiency * HourlyData[(d-1)*48 + t, :Market3Price] for d in days, t in time_hours) -
   sum(charging_power_m3[d,t] / charging_efficiency * HourlyData[(d-1)*48 + t, :Market3Price] for d in days, t in time_hours)
)

# Constraints
## Set the total energy throughput for the whole optimisation horizon 
@constraint(model, cumulative_energy_throughput <= 2 * optimisation_horizon_years * cycles_per_year * max_storage_capacity)

## Set the total energy throughput less than the sum of the total battery charging energy
@constraint(model, cumulative_energy_throughput == 
    sum((charging_power[d, t] * settlement_period_duration) for d in days, t in time_hours) + 
    sum((discharging_power[d, t] * settlement_period_duration) for d in days, t in time_hours))

## Calculate SOE after charging/discharging
@constraint(model, [d in days, t in 1:length(time_hours)-1],  
    storage_level[d, t+1] == storage_level[d, t] + settlement_period_duration * charging_power[d, t] - settlement_period_duration * discharging_power[d, t]
)

## Link the last period of one day to the first period of the next day
@constraint(model, [d in 1:length(days)-1],  
    storage_level[d+1, 1] == storage_level[d, 48] + settlement_period_duration * charging_power[d, 48] - settlement_period_duration * discharging_power[d, 48]
)

## Set the maximum and the minimum storage capacity
@constraint(model, [d in days, t in time_hours], storage_level[d, t] <= max_storage_capacity_limit)
@constraint(model, [d in days, t in time_hours], storage_level[d, t] >= min_storage_capacity_limit)

## Set the initial and the final battery energy level
@constraint(model, [d in days, t in time_hours], storage_level[1, 1] == initial_storage)
@constraint(model, [d in days, t in time_hours], storage_level[end, end] == final_storage)

## Set charging and discharging rate 
@constraint(model, [d in days, t in time_hours], charging_power_m1[d, t] + charging_power_m2[d, t] + charging_power_m3[d,t] == charging_power[d, t])
@constraint(model, [d in days, t in time_hours], discharging_power_m1[d, t] + discharging_power_m2[d, t] + discharging_power_m3[d,t] == discharging_power[d, t])

## Set the maximum charging rate 
@constraint(model, [d in days, t in time_hours], charging_power[d, t] <= max_charge_rate * charge_decision[d, t])
@constraint(model, [d in days, t in time_hours], discharging_power[d, t] <= max_discharge_rate * discharge_decision[d, t])

## Set the maximum charging rate 
@constraint(model, [d in days, t in time_hours], charging_power[d, t] >= min_charge_rate * charge_decision[d, t])
@constraint(model, [d in days, t in time_hours], discharging_power[d, t] >= min_discharge_rate * discharge_decision[d, t])

# Avoid charging and discarging in the same settlement period
@constraint(model, [d in days, t in time_hours], charge_decision[d, t] + discharge_decision[d, t] <= 1)

# Ensure consistent charging or discharging decision for Market 3 across each day (48 time steps)
for d in days
    @constraint(model, charging_power_m3[d, :] .== charging_power_m3[d, 1])
    @constraint(model, discharging_power_m3[d, :] .== discharging_power_m3[d, 1])
 end

# Solve the model
optimize!(model)

# Extract results
println("Extracting results.")
results = DataFrame()

for d in days
    market_prices_m1 = [HourlyData[(d-1)*48 + t, :Market1Price] for t in time_hours]
    market_prices_m2 = [HourlyData[(d-1)*48 + t, :Market2Price] for t in time_hours]
    market_prices_m3 = [HourlyData[(d-1)*48 + t, :Market3Price] for t in time_hours]

    charge_results = [value(charge_decision[d, t]) for t in time_hours]
    discharge_results = [value(discharge_decision[d, t]) for t in time_hours]
    charge_power_results = [value(charging_power[d, t]) for t in time_hours]
    discharge_power_results = [value(discharging_power[d, t]) for t in time_hours]

    charge_m1_power_results = [value(charging_power_m1[d, t]) for t in time_hours]
    discharge_m1_power_results = [value(discharging_power_m1[d, t]) for t in time_hours]

    charge_m2_power_results = [value(charging_power_m2[d, t]) for t in time_hours]
    discharge_m2_power_results = [value(discharging_power_m2[d, t]) for t in time_hours]

    charge_m3_power_results = [value(charging_power_m3[d, t]) for t in time_hours]
    discharge_m3_power_results = [value(discharging_power_m3[d, t]) for t in time_hours]

    storage_results = [value(storage_level[d, t]) for t in time_hours]

    profit_result = 0.5 * [(value(discharging_power_m1[d, t]) * discharging_efficiency - value(charging_power_m1[d, t]) / charging_efficiency) * HourlyData[(d-1)*48 + t, :Market1Price] +
                           (value(discharging_power_m2[d, t]) * discharging_efficiency - value(charging_power_m2[d, t]) / charging_efficiency) * HourlyData[(d-1)*48 + t, :Market2Price] +
                           (value(discharging_power_m3[d, t]) * discharging_efficiency - value(charging_power_m3[d, t]) / charging_efficiency) * HourlyData[(d-1)*48 + t, :Market3Price] for t in time_hours]

    ## Calculate cumulative charged and discharged power
    global cumulative_charged_power += sum(value.(charging_power[d, :])) 
    global cumulative_discharged_power += sum(value.(discharging_power[d, :])) 

    ## Calculate total profits
    total_daily_profit = sum(profit_result)
    global cumulative_profit += total_daily_profit

    ## Calculate the number of cycles
    num_cycles = (cumulative_charged_power + cumulative_discharged_power) / (4 * max_storage_capacity)

    ## Create a DataFrame to store the results of the current day
    daily_results_df = DataFrame(
        Day = fill(unique_dates[d], length(time_hours)),
        Hour = time_hours,
        Market_Price_M1 = market_prices_m1,
        Market_Price_M2 = market_prices_m2,
        Market_Price_M3 = market_prices_m3,

        ### Store the results for the dispatch instruction
        Charge_Decision_Market = charge_results,
        DisCharge_Decision_Market = discharge_results,
        Charge_power_Market = charge_power_results,
        DisCharge_power_Market = discharge_power_results,

        ### Store the results for participation in the first market
        Charge_power_Market1 = charge_m1_power_results,
        DisCharge_power_Market1 = discharge_m1_power_results,

        ### Store the results for participation in the second market
        Charge_power_Market2 = charge_m2_power_results,
        DisCharge_power_Market2 = discharge_m2_power_results,

        ### Store the results for participation in the third market
        Charge_power_Market3 = charge_m3_power_results,
        DisCharge_power_Market3 = discharge_m3_power_results,

        Storage_Level = storage_results,
        Profits = profit_result,
        Cumulative_Profits = fill(cumulative_profit, length(time_hours)),
        Num_Cycles = fill(num_cycles, length(time_hours))
    )
    # Append the daily results to the overall results DataFrame
    append!(results_df, daily_results_df)
end

# Save the combined results to a CSV file
CSV.write("BESS_optimisation_results.csv", results_df)
