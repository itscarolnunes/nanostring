# Analise de expressão gênica 
A plataforma Nanostring nCounter® é amplamente utilizada para pesquisa e aplicações clínicas devido à sua capacidade de medir uma ampla gama de níveis de expressão de mRNA[^1]. A téctinica é altamente sensível, precisa e reprodutível[^2]. Com o intuito de contribuir com a comunidade open source vou compartilhar esse trabalho que realizei e me orgulho muito dele, fiquem a vontade para sugestões!! 😊


A análise é divida em 3 estapas principais que estão esquematizados na imagem ao lado

 <img align= "right" alt="workflow"  src="https://live.staticflickr.com/65535/52647870705_37e1868193_z.jpg">

### 🔹Normalização dos dados 

   A principal e a primeira etapa da análise é a normalização  realizada para reduzir erros aleatórios e efeitos de batches. A precisão e confiabilidade dos resultados dependem da normalização adequada dos dados em relação aos genes de referência. O método de normalização por housekeeping é frequentemente escolhido quando se trata de expressão gênica, entretanto é importante salientar que até mesmo os housekeepings podem sofrer variação na expressão[^3].
   
### 🔹Analise de expressão diferencial 

A analise diferencial é realizada após a normalização para comparar a expressão entre dois grupos. O pacote que utilizo nessa análise é o *limma*. Esse pacote foi desenvolvido primeiramente para análise de Micro Array mas ficou muito popular na análise de RNA-seq. 
O *limma* segue um fluxo de trabalho comum. Primeiro a comparação entre os grupos, que é o desing do experimento. Depois o modelo é ajustado e por fim aplicado o método empirico de Bayes 

O out-put da análise inclui:
- **logFC:** O log do fold-change entre os grupos comparados.
- **t:** A estatística-t para avaliar a expressão diferencial.
- **P.Value:** O valor de P para expressão diferencial, esse valor não está ajustado para múltiplos testes. 
- **adj.P.Val:** O valor de P ajustado para testes múltiplos. Existem diferentes métodos de ajuste o default é o Benjamini-Horchberg.

### 🔹Visualização dos dados

E por fim a visualização dos dados, o tipo gráfico mais utilizado nessa etapa é o heatmap, entretanto existem outras opções, como o volcano plot por exemplo. 
A clusterização é o principal componete da visualização do heatmap, o cluster hieraquico é utilizado nessa análise. 

##

Algumas outras etapas são necessárias para melhorar a análise e obter um resultado mais confiável.
Vou disponibilizar o script de todas essas estapas prévias a normalização. 



#### referências:
[^1]: [A new normalization for Nanostring nCounter gene expression data](https://academic.oup.com/nar/article/47/12/6073/5494770?login=false)
[^2]: [Application of NanoString technologies in companion diagnostic development](https://www.tandfonline.com/doi/full/10.1080/14737159.2019.1623672?scroll=top&needAccess=true&role=tab)
[^3]: [An approach for normalization and quality control for NanoString RNA expression data](https://academic.oup.com/bib/article-abstract/22/3/bbaa163/5891144?login=false)
