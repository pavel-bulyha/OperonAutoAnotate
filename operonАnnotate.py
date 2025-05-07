from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import operondbSearch


def get_gene_dict(record):
    """
    Создает словарь, сопоставляющий имя гена (с приведением к нижнему регистру и удалением подчёркиваний)
    с его фичей из GenBank-записи.
    Рассматриваются фичи типа "gene" или "CDS". Сначала ищется квалификатор "gene",
    затем "locus_tag", а потом "note".
    """
    gene_dict = {}
    for feature in record.features:
        # Обратите внимание: feature.type может быть именно "gene" –
        # используем проверку через in или равенство
        if feature.type == "gene" or feature.type == "CDS":
            gene_name = None
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
            elif "locus_tag" in feature.qualifiers:
                gene_name = feature.qualifiers["locus_tag"][0]
            elif "note" in feature.qualifiers:
                gene_name = feature.qualifiers["note"][0]
            if gene_name:
                key = gene_name.lower().replace("_", "")
                if key not in gene_dict:
                    gene_dict[key] = feature
    return gene_dict


def filter_operon_features(features):
    """
    Принимает список кандидатов на аннотацию (SeqFeature объектов с типом "operon"),
    каждый из которых содержит координаты и квалификаторы.

    Для каждого кандидата вычисляем:
      - start, end, длину операна (end - start)
      - флаг has_info: True, если в квалификаторах присутствует непустой текст по ключам "definition" или "name"
      - порядок (order) — индекс кандидата в исходном списке

    Затем:
      1. Опероны с has_info == True остаются без фильтрации.
      2. Для операнов с has_info == False производится фильтрация:
         если операн A полностью содержится в операне B (B.start <= A.start и A.end <= B.end)
         и длина B строго больше длины A, то A исключается.
         Если длина операнов равна, оставляется тот, чей индекс меньше (то есть встречался раньше).

    Возвращает отфильтрованный список SeqFeature.
    """
    features_with_meta = []
    for i, feat in enumerate(features):
        start = int(feat.location.start)
        end = int(feat.location.end)
        length = end - start
        qualifiers = feat.qualifiers
        has_info = False
        # Если хотя бы один из полей "definition" или "name" непустой — ставим флаг
        if "definition" in qualifiers and qualifiers["definition"].strip():
            has_info = True
        if "name" in qualifiers and qualifiers["name"].strip():
            has_info = True
        features_with_meta.append({
            "index": i,
            "start": start,
            "end": end,
            "length": length,
            "has_info": has_info,
            "feature": feat
        })

    # Разбиваем на две группы:
    groupA = [f for f in features_with_meta if f["has_info"]]
    groupB = [f for f in features_with_meta if not f["has_info"]]

    filteredB = []
    for candidate in groupB:
        keep = True
        for other in groupA + groupB:
            if other["index"] == candidate["index"]:
                continue
            # Проверяем, полностью ли candidate содержится в other
            if other["start"] <= candidate["start"] and candidate["end"] <= other["end"]:
                # Если длина other больше (строго) или если равна, но другой появился раньше, то candidate удаляем
                if (other["length"] > candidate["length"]) or (
                        other["length"] == candidate["length"] and other["index"] < candidate["index"]):
                    keep = False
                    break
        if keep:
            filteredB.append(candidate)

    final = groupA + filteredB
    # Сортируем по исходному порядку (индексу), чтобы сохранить приоритет "первого"
    final_sorted = sorted(final, key=lambda x: x["index"])
    return [item["feature"] for item in final_sorted]


def annotate_operons(genbank_file, output_file, include_all=True):
    """
    Читает GenBank‑файл, извлекает название организма, находит операны через функцию fetch_operon_data,
    определяет рамки операна по первому и последнему гену и добавляет новую операнную фичу.

    Рамки операна определяются как:
      start = минимальная координата первого и последнего гена;
      end = максимальная координата первого и последнего гена.
    Учитывается направление (strand) операна.
    Новая фича получает квалификаторы: operon_id, name, definition, first_gene, last_gene и note.
    """
    # Читаем GenBank‑файл (предполагается, одна запись)
    record = SeqIO.read(genbank_file, "genbank")

    organism = record.annotations.get("organism", "").strip()
    if not organism:
        print("Имя организма не найдено в GenBank файле.")
        sys.exit(1)
    print(f"Обнаружен организм: {organism}")

    # Используем функцию поиска операнов (импортированную из operondbSearch)
    df_operons = operondbSearch.fetch_operon_data(organism, include_all)
    print(f"Найдено операнов: {len(df_operons)}")

    if df_operons.empty:
        print("Опероны не найдены, завершаем.")
        sys.exit(0)

    gene_dict = get_gene_dict(record)

    new_operon_features = []
    for idx, row in df_operons.iterrows():
        genes_list = row["Genes"]
        if not genes_list:
            print(f"Пропускаем строку {idx}: пустой список генов.")
            continue
        # Берем первое и последнее значение из списка генов
        first_gene = genes_list[0]
        last_gene = genes_list[-1]
        # Приводим имена к нижнему регистру и удаляем подчёркивания для сопоставления
        first_gene_key = first_gene.lower().replace("_", "")
        last_gene_key = last_gene.lower().replace("_", "")

        if first_gene_key not in gene_dict:
            print(f"Ген {first_gene} не найден в GenBank для операна, строка {idx}.")
            continue
        if last_gene_key not in gene_dict:
            print(f"Ген {last_gene} не найден в GenBank для операна, строка {idx}.")
            continue

        feature_first = gene_dict[first_gene_key]
        feature_last = gene_dict[last_gene_key]

        start = min(int(feature_first.location.start), int(feature_last.location.start))
        end = max(int(feature_first.location.end), int(feature_last.location.end))
        strand = feature_first.location.strand  # Предполагаем, что операн односторонний

        qualifiers = {}
        qualifiers["operon_id"] = row["Operon ID"]
        # Если есть данные в столбце Name/Definition, добавляем их
        if "Name" in row and row["Name"]:
            qualifiers["name"] = row["Name"]
        if "Definition" in row and row["Definition"]:
            qualifiers["definition"] = row["Definition"]
        qualifiers["first_gene"] = first_gene
        qualifiers["last_gene"] = last_gene
        qualifiers["note"] = f"Operon boundaries: {start}-{end}; strand: {strand}"

        new_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type="operon", qualifiers=qualifiers)
        new_operon_features.append(new_feature)

    if not new_operon_features:
        print("Не удалось создать ни одной аннотации для операнов.")
        sys.exit(0)

    # Фильтруем опероны согласно заданной иерархии:
    filtered_features = filter_operon_features(new_operon_features)
    print(f"После фильтрации осталось {len(filtered_features)} операнов.")

    # Добавляем отфильтрованные операнные фичи к существующим записям в GenBank
    record.features.extend(filtered_features)

    SeqIO.write(record, output_file, "genbank")
    print(f"Аннотированная GenBank запись сохранена в: {output_file}")


if __name__ == '__main__':
    input_gb = "Haemophilus_influenzae_Rd_KW20.gb"  # Замените на путь к вашему GenBank файлу
    output_gb = "annotated.gb"  # Файл для результатов аннотации
    annotate_operons(input_gb, output_gb, include_all=True)