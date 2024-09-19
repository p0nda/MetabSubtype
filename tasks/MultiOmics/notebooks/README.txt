## 使用教程说明文档

### 功能描述

本文档描述了几个用于数据处理和可视化的R语言函数，主要包括以下功能：
1. **统计检验筛选**：通过Wilcoxon检验、t检验和Kruskal-Wallis检验筛选重要特征。
2. **P值计算**：使用Wilcoxon检验和Kruskal-Wallis检验计算特征的P值。
3. **绘制热图**：根据筛选出的特征绘制热图，用于数据的可视化分析。

### 函数列表

1. `filter_by_wilcox`
2. `pvalue_by_wilcox`
3. `pvalue_by_kruskal`
4. `filter_by_ttest`
5. `filter_by_kruskal`
6. `filter_by_pvalue`
7. `draw_heatmap`
8. `draw_single_heatmap`

### 函数详解

#### 1. `filter_by_wilcox`

**功能**：使用Wilcoxon检验筛选显著性特征。

**参数**：
- `df.test`：待检验的数据框。
- `feature_cols`：特征列名向量。
- `class_label`：分类标签列名。
- `threshold`：P值阈值，默认值为100。

**返回**：筛选后的数据框，仅包含显著性特征和分类标签列。

**示例**：
```R
selected_data <- filter_by_wilcox(df.test, feature_cols, class_label, threshold=0.05)
```

#### 2. `pvalue_by_wilcox`

**功能**：计算每个特征的Wilcoxon检验P值。

**参数**：
- `df.test`：待检验的数据框。
- `metab_num`：特征数量。
- `class_label`：分类标签列名。

**返回**：包含特征和对应P值的数据框。

**示例**：
```R
pvalues <- pvalue_by_wilcox(df.test, metab_num, class_label)
```

#### 3. `pvalue_by_kruskal`

**功能**：计算每个特征的Kruskal-Wallis检验P值。

**参数**：
- `df.test`：待检验的数据框。
- `metab_num`：特征数量。
- `class_label`：分类标签列名。

**返回**：包含特征和对应P值的数据框。

**示例**：
```R
pvalues <- pvalue_by_kruskal(df.test, metab_num, class_label)
```

#### 4. `filter_by_ttest`

**功能**：使用t检验筛选显著性特征。

**参数**：
- `df.test`：待检验的数据框。
- `feature_cols`：特征列名向量。
- `class_label`：分类标签列名。
- `threshold`：P值阈值，默认值为100。

**返回**：筛选后的数据框，仅包含显著性特征和分类标签列。

**示例**：
```R
selected_data <- filter_by_ttest(df.test, feature_cols, class_label, threshold=0.05)
```

#### 5. `filter_by_kruskal`

**功能**：使用Kruskal-Wallis检验筛选显著性特征。

**参数**：
- `df.test`：待检验的数据框。
- `feature_cols`：特征列名向量。
- `class_label`：分类标签列名。
- `threshold`：P值阈值，默认值为100。

**返回**：筛选后的数据框，仅包含显著性特征和分类标签列。

**示例**：
```R
selected_data <- filter_by_kruskal(df.test, feature_cols, class_label, threshold=0.05)
```

#### 6. `filter_by_pvalue`

**功能**：根据P值筛选显著性特征，自动选择使用Wilcoxon检验或Kruskal-Wallis检验。

**参数**：
- `df.test`：待检验的数据框。
- `feature_cols`：特征列名向量。
- `class_label`：分类标签列名。
- `threshold`：P值阈值，默认值为100。

**返回**：筛选后的数据框，仅包含显著性特征和分类标签列。

**示例**：
```R
selected_data <- filter_by_pvalue(df.test, feature_cols, class_label, threshold=0.05)
```

#### 7. `draw_heatmap`

**功能**：根据筛选出的特征绘制热图。

**参数**：
- `loaddata`：待绘制的数据框。
- `feature_cols`：特征列名向量。
- `class_label`：分类标签列名。
- `ha_col`：热图注释列名。
- `use_row_ha`：是否使用行注释，默认值为FALSE。
- `col_split`：是否分列显示，默认值为TRUE。

**返回**：绘制的热图对象。

**示例**：
```R
heatmap <- draw_heatmap(loaddata, feature_cols, class_label, ha_col, use_row_ha=TRUE, col_split=TRUE)
```

#### 8. `draw_single_heatmap`

**功能**：根据筛选出的特征绘制单一热图。

**参数**：
- `df.raw`：原始数据框。
- `feature_num`：特征数量。
- `class_label`：分类标签列名。
- `pvalue_cutoff`：P值阈值，默认值为0.05。
- `use_row_ha`：是否使用行注释，默认值为FALSE。

**返回**：绘制的热图对象。

**示例**：
```R
single_heatmap <- draw_single_heatmap(df.raw, feature_num, class_label, pvalue_cutoff=0.05, use_row_ha=TRUE)
```

### 总结

上述函数可以帮助用户对数据进行统计检验筛选，并可视化显著性特征，适用于多种数据分析场景。用户可以根据需要灵活使用这些函数，以获得更具解释性和可视化的数据结果。