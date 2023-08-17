import pandas as pd
import numpy as np
import dask.array as da

def get_hourly_system_capacity_wo_storage(num_iters, for_matrix, nameplate, hourly_capacity):

    outage_status = da.random.random_sample((num_iters, 8760, nameplate.size),
                                        chunks=(150, 150, nameplate.size)) \
                    > for_matrix

    hourly_system_capacity_wo_storage = da.multiply(outage_status,hourly_capacity).sum(axis=2)

    hourly_system_capacity_wo_storage = hourly_system_capacity_wo_storage.persist().rechunk((50,8760))
    return hourly_system_capacity_wo_storage

def get_hourly_storage_capacity(num_iters, hourly_system_capacity_wo_storage, demand, storage):

    reset_storage(storage)

#    hourly_storage_contribution = da.zeros((num_iters,8760),chunks=(100,8760))
    if(num_iters%50!=0):
        chunk_0 = num_iters
    else:
        chunk_0 = 50

    hourly_storage_contribution = da.map_blocks(get_hourly_storage_contribution,hourly_system_capacity_wo_storage,
                                                demand, storage, chunk_0, chunks=(chunk_0,8760),
                                                dtype=hourly_system_capacity_wo_storage.dtype)

    return hourly_storage_contribution

def get_lolp(num_iters, hourly_system_capacity, demand):

    hourly_system_capacity = hourly_system_capacity.compute()

    lolp = ((hourly_system_capacity < demand).sum(axis=0) / num_iters)

    return lolp,hourly_system_capacity

def calc_lolp_full(num_iters, for_matrix, nameplate, hourly_capacity,
             demand, storage, return_hourly_capacity=False):

    hourly_system_capacity_wo_storage = get_hourly_system_capacity_wo_storage(num_iters,for_matrix,nameplate,
                                                                              hourly_capacity)


    if not(storage["num units"] == 0):
        hourly_storage_capacity = get_hourly_storage_capacity(num_iters,
                                                              hourly_system_capacity_wo_storage,
                                                              demand, storage)

        hourly_system_capacity = hourly_system_capacity_wo_storage + hourly_storage_capacity

    else:
        hourly_system_capacity = hourly_system_capacity_wo_storage


    lolp,hourly_system_capacity = get_lolp(num_iters,hourly_system_capacity,demand)

    if return_hourly_capacity:
        return lolp,hourly_system_capacity

    return lolp

def initialize_storage_from_df(storage_df,
                               power_cap_colname,energy_cap_colname,
                               forced_outage_rate,dispatch_strategy,
                               round_trip_efficiency):
# Fill data structure
    storage = dict()
    storage["dispatch strategy"] = dispatch_strategy
    storage["num units"] = len(storage_df)
    storage["max charge rate"] = storage_df[power_cap_colname].values
    storage["max discharge rate"] = storage_df[power_cap_colname].values
    storage["max energy"] = storage_df[energy_cap_colname].values

    # Parametrized (for now)
    storage["roundtrip efficiency"] = np.ones(storage["num units"]) * round_trip_efficiency
    storage["one way efficiency"] = storage["roundtrip efficiency"] ** .5

    # Hourly Tracking (storage starts empty)
    storage["power"] = np.zeros(storage["num units"])
    storage["extractable energy"] = np.ones(storage["num units"])*storage["max energy"]
    storage["energy"] = storage["extractable energy"] / storage["one way efficiency"]
    storage["time to discharge"] = storage["extractable energy"] / storage["max discharge rate"]
    storage["efor"] = forced_outage_rate
    storage["full"] = True

    return storage

def make_storage(include_storage, energy_capacity, charge_rate, discharge_rate, round_trip_efficiency,
                 efor, dispatch_strategy):
# make a storage unit in properly formatted dictionary

    if include_storage == False or energy_capacity == 0:
        storage = dict()
        storage["num units"] = 0
        return storage

    storage = dict()

    storage["dispatch strategy"] = dispatch_strategy
    storage["num units"] = 1
    storage["max charge rate"] = np.array(charge_rate)
    storage["max discharge rate"] = np.array(discharge_rate)
    storage["max energy"] = np.array(energy_capacity)
    storage["roundtrip efficiency"] = np.ones(storage["num units"]) * round_trip_efficiency
    storage["one way efficiency"] = storage["roundtrip efficiency"] ** .5
    storage["power"] = np.zeros(storage["num units"])
    storage["extractable energy"] = np.ones(storage["num units"])*storage["max energy"]
    storage["energy"] = storage["extractable energy"] / storage["one way efficiency"]
    storage["time to discharge"] = storage["extractable energy"] / storage["max discharge rate"]
    storage["efor"] = efor
    storage["full"] = True

    return storage

def append_storage(fleet_storage, additional_storage):
    """ Combine two storage dictionaries

    ...

    Args:
        -----------
        `fleet_storage` (dict): dictionary of storage units
        `additional_storage` (dict): dictionary of storage units
    """
    #edge cases
    if fleet_storage["num units"] == 0:
        return additional_storage

    if additional_storage["num units"] == 0:
        return fleet_storage

    #combine regular attribute(s)
    fleet_storage["num units"] += additional_storage["num units"]

    #concatenate all array objects
    for key in fleet_storage:
        if isinstance(fleet_storage[key], np.ndarray):
            fleet_storage[key] = np.append(fleet_storage[key],additional_storage[key])

    return fleet_storage

def reset_storage(storage):
# All units begin full at the beginning of the year
    #for simulation begin empty
    if storage["num units"] == 0:
        return

    storage["power"] = np.zeros(storage["num units"])
    storage["extractable energy"] = np.ones(storage["num units"])*storage["max energy"] # storage begins full
    storage["energy"] = storage["extractable energy"] / storage["one way efficiency"]
    storage["time to discharge"] = storage["extractable energy"] / storage["max discharge rate"]
    storage["full"] = True

    return

def get_hourly_storage_contribution(hourly_capacity, hourly_load, storage, chunk_0):
    """ Find the hourly capacity matrix for a set of storage units a given number of iterations.

        ...

        Args:
        ----------
        `num_iterations` (int): Number of capacity curves to sample for MCS.

        `hourly_capacity` (ndarray): array of hourly capacities for the desired number of iterations

        `hourly_load` (ndarray): vector of hourly load

        `storage` (dict): dictionary of storage units

        `renewable profile` (ndarray): vector of hourly renewable profiles
    """
    hourly_storage_contribution = np.zeros((chunk_0,8760))

    # edge case
    if storage["num units"] == 0:
        return 0

    # dispatch in every iteration according to outages and available capacity
    for i in range(chunk_0):

        if storage["dispatch strategy"] == "reliability":
            reliability_strategy(hourly_capacity[i,:],hourly_load,storage,hourly_storage_contribution[i,:])

        elif storage["dispatch strategy"] == "arbitrage":
            renewable_profile = None
            net_load = hourly_load - renewable_profile
            arbitrage_strategy(net_load,storage,hourly_storage_contribution[i,:])

        else:
            error_message = "Invalid dispatch strategy: \""+storage["dispatch strategy"]+"\""
            raise RuntimeWarning(error_message)
        reset_storage(storage)

    return hourly_storage_contribution

def arbitrage_strategy(net_load,storage,hourly_storage_contribution):
# BAD IMPLEMENTATION : emulate arbitrage for storage dispatch policy
    for day in range(365):
        start = day*24
        end = (day+1)*24
        arbitrage_dispatch(net_load[start:end],storage,hourly_storage_contribution[start:end])

    return

def arbitrage_dispatch(net_load, storage, hourly_storage_contribution):
# BAD IMPLEMENTATION : shoddy arbitrage emulation with percentile-based peak shaving
    discharge_percentile = 75
    charge_percentile = 25

    discharge_threshold = np.percentile(net_load, discharge_percentile)
    charge_threshold = np.percentile(net_load, charge_percentile)

    for hour in range(24):
        if net_load[hour] > discharge_threshold:
            load_difference = net_load[hour] - discharge_threshold
            hourly_storage_contribution[hour] = discharge_storage(load_difference, storage)
        elif net_load[hour] < charge_threshold:
            load_difference = charge_threshold - net_load[hour]
            hourly_storage_contribution[hour] = charge_storage(load_difference, storage)

    return

def reliability_strategy(hourly_capacity,hourly_load,storage,hourly_storage_contribution):
# Charge/Discharge storage to greedily maximize reliability
    simulation_days = np.unique((np.argwhere(hourly_load > hourly_capacity)//24).flatten())
    simulation_days = np.unique(np.minimum(np.maximum(simulation_days,0),364))

    # no risk days
    if len(simulation_days) != 0:
        last_day = simulation_days[-1]

    #choose strategy on a daily basis
    for i in range(simulation_days.size):
        day = simulation_days[i]
        if day == last_day:
            next_risk_day = day + 1
        else:
            next_risk_day = simulation_days[i+1]
        # Simulate risk days or until storage is charged
        storage["full"] = False
        while day != next_risk_day  and storage["full"] == False:
            start = (day)*24
            end = (day+1)*24
            reliability_dispatch(hourly_storage_contribution[start:end], hourly_load[start:end],
                                hourly_capacity[start:end], storage)
            day += 1

    return

def reliability_dispatch(hourly_storage_contribution, hourly_load, hourly_capacity, storage):
# Disharge when load exceeds capacity, charge otherwise
    floating_point_buffer = 1e-6 # 1 W buffer to avoid false loss-of-loads

    for hour in range(24):
        # discharge if load is not met
        if hourly_load[hour] > hourly_capacity[hour]:
            unmet_load = hourly_load[hour] - hourly_capacity[hour] + floating_point_buffer
            # discharge storage
            hourly_storage_contribution[hour] = discharge_storage(unmet_load, storage)
        # charge if surplus
        else:
            additional_capacity = hourly_capacity[hour] - hourly_load[hour] - floating_point_buffer
            # charge storage
            hourly_storage_contribution[hour] = charge_storage(additional_capacity, storage)

    return

def discharge_storage(unmet_load, storage):
# Discharge according to optimal policy proposed by Evans et.al.
    P_r = unmet_load
    p = storage["max discharge rate"]
    x = storage["time to discharge"]
    u = storage["power"]

    # sort unique "time to discharges", discharge HIGHEST first
    y = np.flip(np.unique(np.concatenate((x,np.maximum(x-1,0)))))
    E_max = 0
    i = 0

    #Find devices partaking in discharge
    while E_max < P_r and i < len(y)-1:
        i += 1
        E_min = E_max
        E_max = np.sum(p*np.maximum(np.minimum(x-y[i],1),0))

    #Load exactly met or unmeetable load
    if E_max <= P_r:
        z = y[i]

    #interpolate
    else:
        z = y[i-1] + (P_r - E_min)/(E_max - E_min)*(y[i]-y[i-1])


    u = p*np.maximum(np.minimum(x-z,1),0) * (np.random.random_sample(x.shape)>storage["efor"])

    if storage["efor"] > 0:
        u *= np.random.random_sample(x.shape)>storage["efor"]

    #update storage
    storage["power"] = u
    update_storage(storage,"discharge")

    return np.sum(storage["power"])

def charge_storage(additional_capacity, storage):
# Charge according to policy proposed by Evans et.al.
    P_r = -additional_capacity
    p_c = -storage["max charge rate"]
    p_d = storage["max discharge rate"]
    x_max = storage["max energy"] / storage["max discharge rate"]
    x = storage["time to discharge"]

    u = storage["power"]

    n = storage["roundtrip efficiency"]

    # sort unique "time to discharges", charge LOWEST first
    z_max = np.minimum(x-n*p_c/p_d, x_max)
    y = np.unique(np.concatenate((x,z_max)))

    E_max = 0
    E_min = 0
    i = 0

    # Find devices partaking in charge
    while E_max < -P_r and i < len(y)-1:
        i += 1
        E_min = E_max
        E_max = np.sum(p_d/n*np.maximum(np.minimum(y[i],z_max)-x,0))

    if E_max <= -P_r:
        z = y[i]
    #interpolate
    else:
        if E_max == 0 and E_min == 0:
            return 0
        else:
            z = y[i-1] + (-P_r - E_min)/(E_max - E_min)*(y[i]-y[i-1])


    u = -1*p_d/n*np.maximum(np.minimum(z,z_max)-x,0)

    if storage["efor"] > 0:
        u *= np.random.random_sample(x.shape)>storage["efor"]

    #update storage
    storage["power"] = u
    update_storage(storage, "charge")

    return np.sum(storage["power"])

def update_storage(storage, status):
# Update charge status of all storage devices

    if status == "discharge":
        storage["extractable energy"] = storage["extractable energy"] - storage["power"]
        storage["energy"] = storage["extractable energy"] / storage["one way efficiency"]
    if status == "charge":
        storage["energy"] = storage["energy"] - storage["power"] * storage["one way efficiency"]
        storage["extractable energy"] = storage["energy"] * storage["one way efficiency"]

    # set storage state
    storage["full"] = np.sum(storage["extractable energy"]) == np.sum(storage["max energy"])

    storage["time to discharge"] = np.divide(storage["extractable energy"],storage["max discharge rate"])

    return

def get_benchmark_fors(benchmark_FORs_file):
    """ Load in benchmark fors for temperature increments of 5 celsius from -15 to 35 for 6 different types of technology

    ...

    Args:
    ----------
    `benchmark_FORs_file` (str): file path to excel file containing temperature dependent fors for different generator technologies
    """

    tech_categories = ["Temperature","HD","CC","CT","DS","HD","NU","ST","Other"]
    forData = pd.read_excel(benchmark_FORs_file, engine='openpyxl')
    benchmark_fors_tech = dict()
    for tech in tech_categories:
        benchmark_fors_tech[tech] = forData[tech].values
    return benchmark_fors_tech
