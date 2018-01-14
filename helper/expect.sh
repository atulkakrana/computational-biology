#!/usr/bin/expect

spawn ssh kakrana@dashi.ddpsc.org
expect "password"
send "dl4sbl1645\r"
interact
