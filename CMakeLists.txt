cmake_minimum_required(VERSION 4.3)
project(seuif97 
    VERSION 0.1.0 
    DESCRIPTION "SEUIF97 - Steam and Water Properties Library"
    LANGUAGES C
)


if(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Set to either \"Release\" or \"Debug\"")
else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Set to either \"Release\" or \"Debug\"")
endif()

# ===== Настройки C стандарта =====
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# ===== Определение платформы =====
if(ANDROID)
    message(STATUS "Building for Android (API level: ${ANDROID_PLATFORM})")
    set(PLATFORM_ANDROID TRUE)
elseif(LINUX)
    message(STATUS "Building for Linux")
    set(PLATFORM_LINUX TRUE)
elseif(WIN32)
    message(STATUS "Building for Windows")
    set(PLATFORM_WINDOWS TRUE)
endif()

# ===== Сбор исходников =====
set(SOURCE_DIRS
    src/common
    src/algo
    src/r1
    src/r2
    src/r3
    src/r4
    src/r5
)

foreach(dir ${SOURCE_DIRS})
    file(GLOB_RECURSE SRC_${dir} CONFIGURE_DEPENDS 
        ${dir}/*.c 
        ${dir}/*.h
    )
    list(APPEND ALL_SOURCES ${SRC_${dir}})
endforeach()

# Фильтрация только .c файлов для компиляции
foreach(src ${ALL_SOURCES})
    if(${src} MATCHES ".*\\.c$")
        list(APPEND LIB_SOURCES ${src})
    endif()
endforeach()

# ===== Настройки библиотеки =====
add_library(seuif97 SHARED)

target_sources(seuif97 PRIVATE ${LIB_SOURCES})

target_include_directories(seuif97 
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/common>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${SOURCE_DIRS}
)

target_compile_definitions(seuif97 
    PRIVATE LIB_EXPORT 
    INTERFACE LIB_IMPORT
)

# ===== Платформозависимые настройки =====
if(PLATFORM_WINDOWS)
    if(MSVC)
        target_compile_options(seuif97 PRIVATE 
            /W3
        )
    elseif(MINGW)
    endif()
elseif(PLATFORM_LINUX OR PLATFORM_ANDROID)
    
    target_compile_options(seuif97 PRIVATE 
        -Wall
        -Wextra
        -fPIC
        -fvisibility=hidden
    )
    # Оптимизация размера для Android
    if(PLATFORM_ANDROID)
        target_compile_options(seuif97 PRIVATE 
            -fdata-sections
            -ffunction-sections
        )
        target_link_options(seuif97 PRIVATE
        "LINKER:--gc-sections"
        "LINKER:--strip-all"
        )
    endif()
endif()

# ===== Демо-приложение =====
add_executable(demo EXCLUDE_FROM_ALL demo/demo-c/demo.c)

target_link_libraries(demo PRIVATE seuif97)

# ===== Установка =====
include(GNUInstallDirs)

install(TARGETS seuif97
    EXPORT seuif97-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Экспорт .lib для Windows
if(PLATFORM_WINDOWS AND MSVC)
    install(FILES 
        ${CMAKE_CURRENT_BINARY_DIR}/seuif97.lib
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        OPTIONAL
    )
endif()

# ===== Android NDK специфичные настройки =====
if(PLATFORM_ANDROID)
    # Настройка ABI
    message(STATUS "Android ABI: ${ANDROID_ABI}")
    message(STATUS "Android STL: ${ANDROID_STL}")
    
    # Минимальная версия API
    if(NOT ANDROID_PLATFORM)
        set(ANDROID_PLATFORM android-21)
    endif()
endif()

# ===== Копирование библиотек =====
add_custom_command(TARGET demo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:seuif97>"
    "$<TARGET_FILE_DIR:demo>"
    COMMENT "Copying library to executable directory"
)

# ===== Дополнительные цели =====

# Установка с проверкой
add_custom_target(install_check
    COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR} --config ${CMAKE_BUILD_TYPE}
    COMMAND ${CMAKE_COMMAND} -E echo "Installation complete. Checking files..."
    COMMAND ${CMAKE_COMMAND} -E echo "Library: $<TARGET_FILE:seuif97>"
)

# ===== Информация о сборке =====
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")