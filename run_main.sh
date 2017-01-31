#!/bin/sh
nohup python3 main.py 1 group_19/ &
nohup python3 main.py 1 group_20/ &
nohup python3 main.py 1 group_21/ &
nohup python3 main.py 1 group_22/ &
nohup python3 main.py 1 group_23/ &
nohup python3 main.py 1 group_24/ &

wait

nohup python3 main.py 1 group_19/ &
nohup python3 main.py 1 group_20/ &
nohup python3 main.py 1 group_21/ &
nohup python3 main.py 1 group_22/ &
nohup python3 main.py 1 group_23/ &
nohup python3 main.py 1 group_24/ &
wait