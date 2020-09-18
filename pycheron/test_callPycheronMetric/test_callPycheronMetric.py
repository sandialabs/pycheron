#####################################################################################
# Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.
#####################################################################################
# NOTICE:
# For five (5) years from 10/21/2019 the United States Government is granted for
# itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and perform publicly and
# display publicly, by or on behalf of the Government. There is provision for the
# possible extension of the term of this license. Subsequent to that period or any
# extension granted, the United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
# data to reproduce, prepare derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so. The specific
# term of the license can be identified by inquiry made to National Technology and
# Engineering Solutions of Sandia, LLC or DOE. NEITHER THE UNITED STATES GOVERNMENT,
# NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING
# SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR
# IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
# USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. Any licensee of this software
# has the obligation and responsibility to abide by the applicable export control laws,
# regulations, and general prohibitions relating to the export of technical data.
# Failure to obtain an export control license or other authority from the Government
# may result in criminal liability under U.S. laws.
# (End of Notice)
####################################################################################

import yaml
import pytest
from pycheron.callPycheronMetric import callPycheron
from pycheron.db.sqllite_db import Database


def load_config_file(config_file):
    with open(config_file, "r") as ymlfile:
        cfg = yaml.load(ymlfile)
    return cfg


def ensure_config_values(config):
    for key, val in list(config.items()):
        if val == "None" and key != "session":
            config[key] = None


def get_db(db_name, session):
    return Database(db_name, session_name=session)


def basic_stats_test(data_base):
    basic_stats = data_base.get_metric("basicStatsMetric")
    assert basic_stats is not None


@pytest.mark.parametrize(
    "config_file",
    [
        "test_callPycheronMetric/pycheronConfigTemplateDir.yaml",
        "test_callPycheronMetric/pycheronConfigTemplateDirNotByDay.yaml",
        "test_callPycheronMetric/pycheronConfigTemplateCalcIndNoPlot.yaml",
        "test_callPycheronMetric/pycheronConfigTemplateCalcIndPlot.yaml",
        "test_callPycheronMetric/pycheronConfigTemplateCalcAllPlot.yaml",
        "test_callPycheronMetric/pycheronConfigTemplateCalcAllNoPlot.yaml",
        "test_callPycheronMetric/pycheronConfigTemplate.yaml",
        "test_callPycheronMetric/pycheronConfigTemplateWfdisc.yaml",
    ],
)
def test_callPycheronMetric(config_file):
    cfg = load_config_file(config_file)
    ensure_config_values(cfg)

    callPycheron(
        output_dir=cfg["output_dir"],
        data=cfg["data"],
        datatype=cfg["datatype"],
        calcAll=cfg["calcAll"],
        calcPsds=cfg["calcPsds"],
        calcBasic=cfg["calcBasic"],
        calcCorr=cfg["calcCorr"],
        calcCrossCorr=cfg["calcCrossCorr"],
        calcGap=cfg["calcGap"],
        calcAmp=cfg["calcAmp"],
        calcSNR=cfg["calcSNR"],
        calcSOH=cfg["calcSOH"],
        calcStalta=cfg["calcStalta"],
        calcDcOffset=cfg["calcDcOffset"],
        calcSpikes=cfg["calcSpikes"],
        calcAllDeadChan=cfg["calcAllDeadChan"],
        calcTransfer=cfg["calcTransfer"],
        calcCal=cfg["calcCal"],
        network=cfg["network"],
        station=cfg["station"],
        byDay=cfg["byDay"],
        startdate=cfg["startdate"],
        enddate=cfg["enddate"],
        jul_start=cfg["jul_start"],
        jul_end=cfg["jul_end"],
        generateMasks=cfg["generateMasks"],
        masksByTime=cfg["masksByTime"],
        rmsThreshold=cfg["rmsThreshold"],
        maxThreshold=cfg["maxThreshold"],
        minThreshold=cfg["minThreshold"],
        medianThreshold=cfg["medianThreshold"],
        meanThreshold=cfg["meanThreshold"],
        varianceThreshold=cfg["varianceThreshold"],
        stdThreshold=cfg["stdThreshold"],
        maxLagSecs=cfg["maxLagSecs"],
        filt=cfg["filt"],
        freqmin=cfg["freqmin"],
        freqmax=cfg["freqmax"],
        corners=cfg["corners"],
        zerophase=cfg["zerophase"],
        maxorder=cfg["maxorder"],
        ba=cfg["ba"],
        freq_passband=cfg["freq_passband"],
        windowSecs=cfg["windowSecs"],
        incrementSecs=cfg["incrementSecs"],
        threshold=cfg["threshold"],
        separateMasks=cfg["separateMasks"],
        completeDay=cfg["completeDay"],
        expLoPeriod=cfg["expLoPeriod"],
        expHiPeriod=cfg["expHiPeriod"],
        linLoPeriod=cfg["linLoPeriod"],
        linHiPeriod=cfg["linHiPeriod"],
        evalresp=cfg["evalresp"],
        dcExpThreshold=cfg["dcExpThreshold"],
        pctBelowNoiseThreshold=cfg["pctBelowNoiseThreshold"],
        pctAboveNoiseThreshold=cfg["pctAboveNoiseThreshold"],
        dcLinThreshold=cfg["dcLinThreshold"],
        num_gaps=cfg["num_gaps"],
        pctBelowNoiseThresholdRESP=cfg["pctBelowNoiseThresholdRESP"],
        pctAboveNoiseThresholdRESP=cfg["pctAboveNoiseThresholdRESP"],
        minRep=cfg["minRep"],
        algorithmSNR=cfg["algorithmSNR"],
        windowSecsSNR=cfg["windowSecsSNR"],
        snrThreshold=cfg["snrThreshold"],
        data_quality=cfg["data_quality"],
        activity=cfg["activity"],
        io_clock=cfg["io_clock"],
        windowSize=cfg["windowSize"],
        thresholdSpikes=cfg["thresholdSpikes"],
        selectivity=cfg["selectivity"],
        fixedThreshold=cfg["fixedThreshold"],
        staSecs=cfg["staSecs"],
        ltaSecs=cfg["ltaSecs"],
        increment=cfg["increment"],
        algorithmSTA=cfg["algorithmSTA"],
        plots=cfg["plots"],
        pdfModel=cfg["pdfModel"],
        per_arr=cfg["per_arr"],
        showNoiseModel=cfg["showNoiseModel"],
        showMaxMin=cfg["showMaxMin"],
        showMode=cfg["showMode"],
        showMean=cfg["showMean"],
        showMedian=cfg["showMedian"],
        showEnvelope=cfg["showEnvelope"],
        envelopeType=cfg["envelopeType"],
        showSingle=cfg["showSingle"],
        singleType=cfg["singleType"],
        min_stations=cfg["min_stations"],
        rank_by=cfg["rank_by"],
        processesPSD=cfg["processesPSD"],
        processesSpikes=cfg["processesSpikes"],
        log=cfg["log"],
        fortran=cfg["fortran"],
        timespan=cfg["timespan"],
        dcADF_win_size=cfg["dcADF_win_size"],
        dcADF_pval_thresh=cfg["dcADF_pval_thresh"],
        dcADF_threshold=cfg["dcADF_threshold"],
        dcADF_use_thresh=cfg["dcADF_use_thresh"],
        dcMean_win_size=cfg["dcMean_win_size"],
        dcMean_thresh=cfg["dcMean_thresh"],
        cal_metric_store=cfg["cal_metric_store"],
        dcExpThresholdHour=cfg["dcExpThresholdHour"],
        dcLinThresholdHour=cfg["dcLinThresholdHour"],
        byHourOn=cfg["byHourOn"],
        database=cfg["database"],
        session=cfg["session"],
        overwrite=cfg["overwrite"],
        to_csv=cfg["to_csv"],
        stationStartAt=cfg["stationStartAt"],
    )

    db_name = cfg["output_dir"] + "/" + cfg["database"]
    db = get_db(db_name, cfg["session"])

    basic_stats_test(db)
