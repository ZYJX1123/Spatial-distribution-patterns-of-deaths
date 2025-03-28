import pandas as pd
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import random
from scipy.stats import mannwhitneyu

# 计算两点之间的欧式距离
def euclidean_distance(p1, p2):
    return np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

# 寻找 Voronoi 邻居并计算距离
def find_voronoi_neighbors(points):
    vor = Voronoi(points)
    neighbors = []
    for ridge_vertices in vor.ridge_vertices:
        if -1 not in ridge_vertices:
            for i, vertex in enumerate(vor.ridge_vertices):
                if np.array_equal(vertex, ridge_vertices):
                    p1_index, p2_index = vor.ridge_points[i]
                    p1 = points[p1_index]
                    p2 = points[p2_index]
                    distance = euclidean_distance(p1, p2)
                    neighbors.append((p1_index, p2_index, distance))
                    break
    return neighbors

# 计算列表的均值
def calc_mean(l):
    return sum(l) / float(len(l))

# 评估真实数据和随机排列数据的均值差异，使用 Mann - Whitney U 检验
def assess_mean(real_distances, permuted_distances_list):
    all_permuted_distances = []
    for distances in permuted_distances_list:
        all_permuted_distances.extend(distances)

    # 执行 Mann - Whitney U 检验
    _, p_value = mannwhitneyu(real_distances, all_permuted_distances)

    real_mean_distance = calc_mean(real_distances)
    scramble_mean_distance = calc_mean(all_permuted_distances)

    permuted_means = [calc_mean(distances) for distances in permuted_distances_list]
    permuted_means.sort()
    n_scramble = len(permuted_distances_list)
    scramble_signficance_95 = permuted_means[int(n_scramble * 5 / 100)]
    scramble_signficance_98 = permuted_means[int(n_scramble * 2 / 100)]

    return p_value, real_mean_distance, scramble_mean_distance, scramble_signficance_95, scramble_signficance_98

# 读取 Excel 文件
file_path = r"D:\cell1\2.csv"  # 替换为你的csv文件路径
df = pd.read_csv(file_path)

# 步骤 1: 选取前 n 个细胞的坐标数据
n = 123 #可根据需要修改 n 的值
first_n_points = df.iloc[:n, 1:3].values
real_neighbors = find_voronoi_neighbors(first_n_points)
real_neighbors_df = pd.DataFrame(real_neighbors, columns=['Cell1', 'Cell2', 'Distance'])
real_neighbors_df.to_csv('real_neighbors_stats.csv', index=False)
real_distances = real_neighbors_df['Distance'].tolist()

# 步骤 2: 随机选取 n 个细胞的坐标数据，循环 1000 次
n_scramble = 1000
all_permuted_distances = []
permuted_distances_list = []
for _ in range(n_scramble):
    random_points = df.sample(n).iloc[:, 1:3].values
    permuted_neighbors = find_voronoi_neighbors(random_points)
    distances = [neighbor[2] for neighbor in permuted_neighbors]
    all_permuted_distances.extend(distances)
    permuted_distances_list.append(distances)

permuted_neighbors_df = pd.DataFrame({'Distance': all_permuted_distances})
permuted_neighbors_df.to_csv('permuted_neighbors_stats.csv', index=False)

# 步骤 3: 数据统计分析
m = 20  # 可根据需要修改 m 的值
min_distance = permuted_neighbors_df['Distance'].min()
max_distance = permuted_neighbors_df['Distance'].max()
bins = np.linspace(min_distance, max_distance, m + 1)
hist, _ = np.histogram(permuted_neighbors_df['Distance'], bins=bins)
total_count = len(permuted_neighbors_df['Distance'])
proportions = hist / total_count
stats_analyzed_df = pd.DataFrame({
    'Bin_Min': bins[:-1],
    'Bin_Max': bins[1:],
    'Count': hist,
    'Proportion': proportions
})
stats_analyzed_df.to_csv('permuted_neighbors_stats_analyzed.csv', index=False)

# 步骤 4: 比较真实数据和随机排列数据的统计差异
p_value, real_mean, scramble_mean, scramble_signficance_95, scramble_signficance_98 = assess_mean(real_distances, permuted_distances_list)

# 输出统计结果
print(f"p-value: {p_value}")
print(f"Real mean distance: {real_mean}")
print(f"Scrambled mean distance: {scramble_mean}")
print(f"Scramble significance 95%: {scramble_signficance_95}")
print(f"Scramble significance 98%: {scramble_signficance_98}")

# 将统计结果保存到 CSV 文件
stats_comparison_df = pd.DataFrame({
    'p_value': [p_value],
    'real_mean_distance': [real_mean],
    'scramble_mean_distance': [scramble_mean],
    'scramble_significance_95': [scramble_signficance_95],
    'scramble_significance_98': [scramble_signficance_98]
})
stats_comparison_df.to_csv('stats_comparison.csv', index=False)