git init


:: git config --global user.email "andreysmorodov@gmail.com"
:: git config --global user.name "Smorodov"
@echo off
echo Инициализация Git репозитория...
git init

echo Добавление всех файлов...
git add .

echo Создание первого коммита...
git commit -m "Initial commit"

echo Установка ветки main...
git branch -M main

echo Добавление удалённого репозитория origin...
git remote add origin https://github.com/Smorodov/bell_design.git

echo Готово! Запустите git-remote-add-origin.bat если нужно изменить remote
pause
