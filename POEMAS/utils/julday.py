from datetime import datetime, timedelta

def time(jd):

    time = str((datetime(2001, 1, 1) + timedelta(seconds=jd)).time())
    return time

def date(jd):

    date = str((datetime(2001, 1, 1) + timedelta(seconds=jd)).date())
    return date

def msec(sec):

    msec = timedelta(seconds=sec).seconds * 1000
    return msec