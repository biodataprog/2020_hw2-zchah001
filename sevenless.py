#Print out values from 0 to 99 except those  perfectly divisible by 7
#!/usr/bin/env python3

Start = 0
End   = 99
divisor=7
print("Print out  numbers from",Start,"to",End, " not divisible by",divisor)
for Sevenless  in range(0, 99):
     if Sevenless  % 7 != 0:
                print(Sevenless)
done

