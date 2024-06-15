import sys
import numpy as np
import math
import matplotlib.pyplot as plt
def build_trinomial_tree(T, N, r0, theta, kappa, sigma, Wiener, updownFactor): 
    dt = T / N
# check if sigma value is valid
    if 2 * kappa * theta < sigma ** 2:
        print(f"sigma:{sigma} kappa:{kappa} theta:{theta} Please decrease sigma:{sigma}")
        sys.exit()
    # Probabilities of going up and down:
    p_up = ((math.exp(r0 * dt / 2) - math.exp(-sigma * math.sqrt(dt / 2))) / ( math.exp(sigma * math.sqrt(dt / 2)) - math.exp(-sigma * math.sqrt(dt / 2)))) ** 2
    p_down = ((math.exp(sigma * math.sqrt(dt / 2)) - math.exp(r0 * dt / 2)) / ( math.exp(sigma * math.sqrt(dt / 2)) - math.exp(-sigma * math.sqrt(dt / 2)))) ** 2
    p_middle = 1 - p_up - p_down
    if debugging:
        print(f"p_up:{p_up} p_down:{p_down} p_middle:{p_middle}")
    if p_up >= 1 or p_down >= 1 or p_up + p_down >= 1:
        print(f"p_up:{p_up} p_down:{p_down}! Please increase sigma:{sigma}") 
        sys.exit()
    r_up = math.exp(+sigma * math.sqrt(2 * dt)) * updownFactor 
    r_down = math.exp(-sigma * math.sqrt(2 * dt)) * updownFactor 
    if debugging:
        print(f"p_up:{p_up} p_down:{p_down} p_middle:{p_middle}")
    if p_up >= 1 or p_down >= 1 or p_up + p_down >= 1:
        print(f"p_up:{p_up} p_down:{p_down}! Please increase sigma:{sigma}") 
        sys.exit()
    if debugging:
        print(f"	r_up:{r_up} r_down:{r_down}")
    r_tree = np.zeros((N + 1, 2 * N + 1))
    r_tree[0, N] = r0
    for i in range(N):
        if not debugging:
            print(f"Building trinomial tree on step {i}")
        for j in range(N - i, N + i + 1, 1):
            random_factor = 0
            if Wiener:
                random_factor = sigma * np.sqrt(r_tree[i, j]) * np.random.normal(loc=0, scale=np.sqrt(dt)) 
                if debugging:
                    print(f"random_factor:{random_factor}")
            delta_r = kappa * (theta - r_tree[i, j]) * dt + random_factor 
            if debugging:
                print(f"i:{i} j:{j} delta_r:{delta_r}")

            if i == 0:
                r_tree[i + 1, j] = r_tree[i, j] + delta_r
                r_tree[i + 1, j + 1] = r_tree[i, j] + r_up + delta_r 
                r_tree[i + 1, j - 1] = r_tree[i, j] - r_down + delta_r 
            else:
                if j == N + i: # right side of the tree
                    r_tree[i + 1, j + 1] = r_tree[i, j] + r_up + delta_r
                    r_tree[i + 1, j] = r_tree[i, j] * p_middle /(p_middle + p_up) + (r_tree[i, j - 1]+r_up) * p_up / (p_middle + p_up) + delta_r
                elif j == N - i: # left side of the tree
                    r_tree[i + 1, j - 1] = r_tree[i, j] - r_down + delta_r
                    r_tree[i + 1, j] = r_tree[i, j] * p_middle /(p_middle + p_down) + (r_tree[i, j + 1]-r_down) * p_down / (p_middle + p_down) + delta_r
                else:
                    r_tree[i + 1, j] = r_tree[i, j] * p_middle + (r_tree[i, j - 1]+r_up) * p_up + (r_tree[i, j + 1]-r_down) * p_down + delta_r

    return r_tree
def draw_tree(tree, label):
    levels = tree.shape[0]
    if levels > 6:
        label = 0
    for i in range(N, -1, -1):
        for j in range(N - i, N + i + 1, 1):
            plt.plot(j, i, 'o', markersize=10, color='orange')
            if label:
                plt.text(j, i, f'{tree[i, j]:.4f}', ha='center', va='center', fontsize=8, color='black')
            if i < N:
                plt.plot([j, j + 1], [i, i + 1], color='orange') # Connect to the upper        -right node
                plt.plot([j, j], [i, i + 1], color='orange') # Connect to the upper        node
                plt.plot([j, j - 1], [i, i + 1], color='orange') # Connect to the upper        -left node
    plt.xlabel('IR change (up/down)')
    plt.ylabel('Time Step')
    plt.title('Trinomial Tree for Cox Ingersoll Ross Model')
    plt.show()
def draw_average(tree, N):
    avg_rates = []
    levels = tree.shape[0]
    for i in range(levels):
        average = 0
        count = 0
        if debugging:
            print(f"average rate for i:{i}")
        for j in range(N - i, N + i + 1, 1):
            if debugging:
            	print(f"	j:{j} tree[i, j]:{tree[i, j]}")
            average += tree[i, j]
            count += 1
        avg_rates.append(average / count)
# Plot average interest rate vs time step
    plt.plot(avg_rates, range(len(avg_rates)), marker='o', color='blue') 
    plt.xlabel('Average Interest Rate')
    plt.ylabel('Time Step')
    plt.title('Average Interest Rate vs Time Step')
    if levels < 20:
        for i, txt in enumerate(avg_rates):
            plt.text(txt, i, f"{txt:.3f}", ha='left', va='center', fontsize=8)
    else:
        plt.text(avg_rates[levels-1], levels-1, f"{avg_rates[levels-1]:.3f}", ha='left', va='center', fontsize=8)
    plt.show()
# parameters 
T = 1
N = 50
r0 = 0.05
theta = 0.06
kappa = 2.7 
sigma = 0.3
Wiener = 1
label = 1
updownFactor = 0.001
debugging = 0
tree = build_trinomial_tree(T, N, r0, theta, kappa, sigma, Wiener, updownFactor) 
draw_tree(tree, label)
draw_average(tree, N)

