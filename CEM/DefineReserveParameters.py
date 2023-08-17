# Define reserve parameters for UC & CE w/ UC constraints models
def defineReserveParameters(stoMkts,stoFTLabels):
    # Regulation eligibility
    regElig = ['Steam', 'Combined Cycle', 'Geothermal'] + (stoFTLabels if 'res' in stoMkts.lower() else [])
    contFlexInelig = ['Wind','Solar','DAC'] #fuel types that are ineligible to provide flex or cont reserves

    # Regulation provision cost as fraction of operating cost
    regCostFrac = 0  # 0.1

    # Requirement parameters - based on WWSIS Phase 2
    regLoadFrac = .01                                   # frac of hourly load in reg up & down
    contLoadFrac = .03                                  # frac of hourly load in contingency
    regErrorPercentile = 40                             # ptl of hourly W&S forecast errors for reg reserves; in WWSIS, 95th ptl of 10-m wind & 5-m solar forecast errors
    flexErrorPercentile = 70                            # ptl of hourly W&S forecast errors used in flex reserves

    # Timeframes
    regReserveMinutes = 5               # reg res must be provided w/in 5 minutes
    flexReserveMinutes = 10             # spin reserves must be provided w/in 10 minutes
    contingencyReserveMinutes = 30      # contingency res must be provided w/in 30 minutes
    minutesPerHour = 60
    rampRateToRegReserveScalar = regReserveMinutes/minutesPerHour               # ramp rate in MW/hr
    rampRateToFlexReserveScalar = flexReserveMinutes/minutesPerHour             # ramp rate in MW/hr
    rampRateToContReserveScalar = contingencyReserveMinutes/minutesPerHour

    return (regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile,
        regElig, contFlexInelig, regCostFrac, rampRateToRegReserveScalar, rampRateToFlexReserveScalar,
        rampRateToContReserveScalar)